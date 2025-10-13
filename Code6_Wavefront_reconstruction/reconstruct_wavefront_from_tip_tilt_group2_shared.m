% Reconstruct piston/tip/tilt (PTT) on the IrisAO DM grid from PFAM
% tip/tilt measured in three rounds, then visualize and save wrapped DM commands
%
% Dependencies (in ./function on MATLAB path):
%   Find6NeighbrIdx, CalcLockedSegShift, weight_display, ZernikeDecomposition

%% -------------------------- Configuration ---------------------------------
close all; clear; clc;

%% ------------------------- Folder / Path Setup ----------------------------
currentFolder = pwd;
[~, folderName] = fileparts(currentFolder);
if ~strcmpi(folderName, 'Code6_Wavefront_reconstruction')
    error('Run this script from "Code6_Wavefront_reconstruction". Current folder: %s', folderName);
end

rootFolder = fileparts(currentFolder);
addpath(fullfile(rootFolder, 'function'));          % helper functions
dmDataDir  = fullfile(rootFolder, 'DM_data');    % DM-related lookups
dataDir    = fullfile(currentFolder, 'data');       % local data
grpDir     = fullfile(dataDir, '2_Group');          % round1/2/3 .mat files


%% --------------------------- Parameters -----------------------------------
Param.level = 7;
Param.SegN  = Param.level*(Param.level+1)*3 + 1;   % L7 → 169
Param.UnreachableSeg = [128, 135, 142, 149, 156, 163];
Param.LockedSeg      = [5, 70, 80, 100, 155];
Param.Outlier_hex    = [169,128,129,134,135,136,141,142,143,148,149,150,155,156,157,162,163,164,5,70,80,100];

DM_period   = 77;        % pupil diameter used for display
SHcam_pitch = 6.5;       % µm
DM_pitch    = 606;       % µm

exclude_defocus = 1;     % whether to remove defocus during PTTD removal

%% ------------------------ DM coordinate mappings / Load GT -----------------------
% Expected content:
%   - DM_geometry.mat: PosRef (spot/segment geometry for display)
%   - IndexTransfer.mat: IndexSeq_hex2jk, IndexCell_jk2hex (index conversions)
load(fullfile(dmDataDir, 'DM_geometry.mat'));  % -> PosRef
Xfer = load(fullfile(dmDataDir, 'IndexTransfer.mat'), 'IndexSeq_hex2jk', 'IndexCell_jk2hex');
Convert_hex2jk = Xfer.IndexSeq_hex2jk;
Convert_jk2hex = Xfer.IndexCell_jk2hex;

% Dummy DM image canvas for weight_display (size must be large enough)
DM_im = zeros(1312, 1312);

% Zernike→PTT dictionaries
Z2PTT = load(fullfile(dmDataDir, 'ZernikeCoeff2PTT2.mat'), 'PTT_mode');        % used later to remove Zernike modes
PTT_mode_base  = Z2PTT.PTT_mode;

% Ground truth (astigmatism correction)
IrisAO_PTT_GT = -0.1*squeeze(PTT_mode_base(4,:,:));

% Initial PTT offsets (zeros if none, load the previous value if this is second or third round)
IrisAO_PTT0_1 = zeros(Param.SegN, 3);
IrisAO_PTT0_2 = zeros(Param.SegN, 3);

%% ------------------------ Load tip/tilt per round -------------------------
% Each file must contain a variable "tip_tilt_saved" of size [N x 2] (xTilt,yTilt in radians/deg? kept as-is)
round1_tip_tilt = load(fullfile(grpDir, 'round1_SIFT_ft.mat'), 'tip_tilt_saved').tip_tilt_saved;
round2_tip_tilt = load(fullfile(grpDir, 'round2_SIFT_ft.mat'), 'tip_tilt_saved').tip_tilt_saved;

% Center by grand mean (remove global bias)
round123_center = mean([round1_tip_tilt; round2_tip_tilt;], 1, 'omitnan');
round1_tip_tilt = round1_tip_tilt - round123_center;
round2_tip_tilt = round2_tip_tilt - round123_center;

%% ------------------------ Target segments per round -----------------------
target_segment1 = [1 20 23 26 29 32 35];
target_segment2 = [92 95 98 101 104 107 110 113 116 119 122 125];

%% -------- Propagate tip/tilt angles to each seed and its six neighbors -----------
tip_tilt_total1 = NaN(Param.SegN, 2);
tip_tilt_total2 = NaN(Param.SegN, 2);

round1_segment = [];
round2_segment = [];

for ii = 1:numel(target_segment1)
    seg = target_segment1(ii);
    [neigh, ~] = Find6NeighbrIdx(seg, Param, Convert_hex2jk, Convert_jk2hex);
    if ~(seg == 1 || seg == 2), neigh(neigh==1) = 163; end
    bundle = [seg; neigh];
    round1_segment = [round1_segment; bundle'];
    tip_tilt_total1(bundle, :) = repmat(round1_tip_tilt(ii,:), 7, 1);
end
for ii = 1:numel(target_segment2)
    seg = target_segment2(ii);
    [neigh, ~] = Find6NeighbrIdx(seg, Param, Convert_hex2jk, Convert_jk2hex);
    if ~(seg == 1 || seg == 2), neigh(neigh==1) = 163; end
    bundle = [seg; neigh];
    round2_segment = [round2_segment; bundle']; 
    tip_tilt_total2(bundle, :) = repmat(round2_tip_tilt(ii,:), 7, 1);
end

used_segment_total = unique([round1_segment; round2_segment;]);

%% ---------------------- Macro-segments PTT mapping bases -----------------------
% Overwrite rows 2 (tip) & 3 (tilt) with your Macro-segment PTT values.
PTT_mode_round1 = PTT_mode_base;
PTT_mode_round2 = PTT_mode_base;

row2_block = [ ...
    0.00204724857001837,-1.20486197704087e-18,0.630552559565678; ...
    0.00204724857001847,-4.53257600886802e-17,0.630552559565683; ...
    0.276378556952489,  9.79897033613380e-16,0.630552559565678; ...
    0.276378556952489, -9.47767380892290e-16,0.630552559565678; ...
    0.00204724857001847,0,                  0.630552559565678; ...
   -0.272284059812451,  9.99232199625893e-16,0.630552559565673; ...
   -0.272284059812451, -1.00336315497575e-15,0.630552559565672];
row3_block = [ ...
    0.00204724857001844,0.630552559565684, 2.26820881735948e-17; ...
    0.317323528352860,  0.630552559565684, 9.63736666622372e-16; ...
    0.157638139891421,  0.630552559565682, 9.88959878367924e-16; ...
   -0.157638139891421,  0.630552559565682,-6.92432182390090e-16; ...
   -0.317323528352860,  0.630552559565682, 0; ...
   -0.157638139891421,  0.630552559565682,-9.03454401488104e-17; ...
    0.157638139891421,  0.630552559565682, 4.12648576278340e-16];

for t = 1:numel(target_segment1)
    PTT_mode_round1(2, round1_segment(t,:), :) = row2_block;
    PTT_mode_round1(3, round1_segment(t,:), :) = row3_block;
end
for t = 1:numel(target_segment2)
    PTT_mode_round2(2, round2_segment(t,:), :) = row2_block;
    PTT_mode_round2(3, round2_segment(t,:), :) = row3_block;
end

%% ------------------- Convert tip/tilt angles to PTT -------------------
round1_PTT = tip_tilt_total1(:,1).*squeeze(PTT_mode_round1(2,:,:)) + ...
             tip_tilt_total1(:,2).*squeeze(PTT_mode_round1(3,:,:));
round2_PTT = tip_tilt_total2(:,1).*squeeze(PTT_mode_round2(2,:,:)) + ...
             tip_tilt_total2(:,2).*squeeze(PTT_mode_round2(3,:,:));

PTT_total = zeros(2, Param.SegN, 3);
PTT_total(1,:,:) = round1_PTT;
PTT_total(2,:,:) = round2_PTT;
PTT_total_mean = squeeze(mean(PTT_total, 1, 'omitnan'));

%% -------------------------- Clean results & inpaint missing segments via 6-neighbor interpolation -----------------------------
flag_nan = find(isnan(PTT_total_mean(:,1)));
Param.Outlier_hex = union(Param.Outlier_hex, flag_nan);
PTT_total_mean(flag_nan,:,:) = 0;

unused_segment_total = setdiff(1:Param.SegN, used_segment_total(:)');
unused_segment_total = unique([unused_segment_total, flag_nan']); 

% Regularize locked/unused segments
PTT_total_mean = CalcLockedSegShift(PTT_total_mean, setdiff(1:Param.SegN, used_segment_total(:)'), Param, Convert_hex2jk, Convert_jk2hex);

% Remove segment-1 offset (piston and TT reference)
Seg_active = 1:Param.SegN;
PTT_total_mean(Seg_active, 1) = PTT_total_mean(Seg_active, 1) - PTT_total_mean(1, 1);
PTT_total_mean(Seg_active, 2) = PTT_total_mean(Seg_active, 2) - PTT_total_mean(1, 2);
PTT_total_mean(Seg_active, 3) = PTT_total_mean(Seg_active, 3) - PTT_total_mean(1, 3);

%% ------------------------ Reconstruct phase on hex grid -------------------
slopes = zeros(Param.SegN, 2);
slopes(:,1) = -2/1000 * PTT_total_mean(:, 3);  % x-slope
slopes(:,2) = -2/1000 * PTT_total_mean(:, 2);  % y-slope

slopes_new = slopes;
flag_nan = find(isnan(slopes_new(:,1)));
Param.Outlier_hex = union(Param.Outlier_hex, flag_nan);
slopes_new(flag_nan,:) = 0;

W        = zeros(Param.SegN, 1);   % current wavefront (um)
W_update = zeros(Param.SegN, 1);   % updated wavefront
tol = 1e-3; iterN = 0; delta_W = []; 
Seg_active = 1:Param.SegN; Seg_active(flag_nan) = [];

disp('Looping for phase reconstruction...');
while true
    iterN = iterN + 1;
    if mod(iterN, 100) == 0, disp(iterN); end

    for p = Seg_active
        [neigh, w] = Find6NeighbrIdx(p, Param, Convert_hex2jk, Convert_jk2hex);

        slope_p = [slopes(p,1), slopes(p,2)];      % [sx, sy]
        s_hex   = slopes(neigh, :);                % [sx_i, sy_i]
        sw      = sum(w);

        term0 = (w.' * W(neigh)) / sw;
        term1 = DM_pitch/4 * ( ...
            slope_p(2) * (2*w(1) - 2*w(4) + w(2) - w(3) - w(5) + w(6)) + ...
            slope_p(1) * sqrt(3) * (w(2) + w(3) - w(5) - w(6)) );
        term2 = DM_pitch/4 * ( ...
            (2*w(1)*s_hex(1,2) - 2*w(4)*s_hex(4,2) + w(2)*s_hex(2,2) - w(3)*s_hex(3,2) - w(5)*s_hex(5,2) + w(6)*s_hex(6,2)) + ...
             sqrt(3)*(w(2)*s_hex(2,1) + w(3)*s_hex(3,1) - w(5)*s_hex(5,1) - w(6)*s_hex(6,1)) );

        W_update(p) = term0 + (term1 + term2) / sw;  % um
    end

    err = norm(W_update - W) / max(norm(W_update), eps);
    W = W_update;
    if err < tol, break; end
end
fprintf('Finish phase reconstruction, # of iterations = %d.\n', iterN);

% Remove global piston offset
W_update = W_update - mean(W_update(:));

%% ----------------------- Build DM commands & display ----------------------
IrisAO_PTT_res = zeros(Param.SegN, 3);
IrisAO_PTT_res(:, 1) = -W_update;            % piston (um)
IrisAO_PTT_res(:, 2) = slopes(:, 2) * 1000;  % tip   (mrad)
IrisAO_PTT_res(:, 3) = slopes(:, 1) * 1000;  % tilt  (mrad)

% Reflective DM → divide by 2
IrisAO_PTT_res = -IrisAO_PTT_res / 2;
IrisAO_PTT     = IrisAO_PTT0_1 + IrisAO_PTT_res;

% Display wavefronts
wavefront1         = -2 * IrisAO_PTT_res(:, 1);                    % measured aberration
wavefront2         = -2 * IrisAO_PTT(:, 1);                         % corrective wavefront
wavefront_residual = -2 * (IrisAO_PTT(:, 1) - IrisAO_PTT_GT(:, 1)); % residual

res_slopes        = IrisAO_PTT_res(:, 3:-1:2) * 2 / 1000 * SHcam_pitch;
collective_slopes = IrisAO_PTT(:,  3:-1:2)    * 2 / 1000 * SHcam_pitch;
residual_slopes   = (IrisAO_PTT(:,3:-1:2) - IrisAO_PTT_GT(:,3:-1:2)) * 2 / 1000 * SHcam_pitch;

wavefront1_d        = weight_display(DM_im, PosRef, wavefront1,        DM_period, Param.UnreachableSeg, res_slopes);
wavefront2_d        = weight_display(DM_im, PosRef, wavefront2,        DM_period, Param.UnreachableSeg, collective_slopes);
wavefrontResidual_d = weight_display(DM_im, PosRef, wavefront_residual,DM_period, Param.UnreachableSeg, residual_slopes);

figure(4); set(gcf,'Position',[10 450 1100 300]);
subplot(1,3,1); imagesc(wavefront1_d        - min(wavefront1_d(:)));        axis image; colorbar; colormap jet;  title('Measured Aberration');
subplot(1,3,2); imagesc(wavefront2_d        - min(wavefront2_d(:)));        axis image; colorbar; colormap jet;  title('Corrective Wavefront');
subplot(1,3,3); imagesc(wavefrontResidual_d - min(wavefrontResidual_d(:))); axis image; colorbar; colormap jet; title('Residual Wavefront');

% Zernike decomposition (mask within pupil)
pupilSize   = round((mean(PosRef(158:161,1)) - mean(PosRef(137:140,1))) / Param.level * (Param.level + 1/6));
centerPos(1)= round((mean(PosRef(117:121,1)) + mean(PosRef([99,101:103],1))) / 2);
centerPos(2)= round((mean(PosRef([146,109,111,152],2)) + mean(PosRef([131,93,127,167],2))) / 2);
figure(4), hold on; ang = 0:0.01:2*pi; xp = pupilSize/2*cos(ang); yp = pupilSize/2*sin(ang);
plot(centerPos(1)+xp, centerPos(2)+yp, 'r'); hold off;
centerPos = [size(wavefront1_d,1) - centerPos(2), size(wavefront1_d,2) - centerPos(1)];
zernikeCoeff1 = ZernikeDecomposition(rot90(fliplr(wavefront1_d),        3), centerPos, pupilSize);
zernikeCoeff2 = ZernikeDecomposition(rot90(fliplr(wavefront2_d),        3), centerPos, pupilSize);
zernikeCoeff3 = ZernikeDecomposition(rot90(fliplr(wavefrontResidual_d), 3), centerPos, pupilSize);

figure(5), set(gcf, 'Position', [10 100 800 300]);
subplot(1,3,1); bar(zernikeCoeff1); title('Measured Zernike');  xlabel('Mode'); ylabel('um rms');
subplot(1,3,2); bar(zernikeCoeff2); title('Corrective Zernike');xlabel('Mode'); ylabel('um rms');
subplot(1,3,3); bar(zernikeCoeff3); title('Residual Zernike'); xlabel('Mode'); ylabel('um rms');

% Residual RMS (masked)IrisAO_PTT_gtIrisAO_PTT_gt
residual_wavefront1_zero = (wavefrontResidual_d - min(wavefrontResidual_d(:)));
mask = poly2mask(centerPos(1)+xp, centerPos(2)+yp, size(residual_wavefront1_zero,1), size(residual_wavefront1_zero,2));
rms_error1        = std(residual_wavefront1_zero(mask)); 
rms_error1_lambda = rms_error1 / 0.515;                  

%% ------------------ Remove PTTD (piston, tip, tilt & defocus) -------------
Param.WeakSeg = [169,129,134,136,141,143,148,150,155,157,162,164];
Param.UnreachableSeg = unique([Param.UnreachableSeg, Param.WeakSeg]); % only for display masks
Param.isRmvPTTD = 1;

% Phase wrapping (round 1) for display comparison
phase_thres = 0.515;  % ±lambda [um]
phase_wrap  = IrisAO_PTT(:, 1);
while any(phase_wrap  >  phase_thres), phase_wrap( phase_wrap >  phase_thres) = phase_wrap( phase_wrap >  phase_thres) - phase_thres; end
while any(phase_wrap  < -phase_thres), phase_wrap( phase_wrap < -phase_thres) = phase_wrap( phase_wrap < -phase_thres) + phase_thres; end

if Param.isRmvPTTD
    % Project-out Zernike modes 1–3 (& 5 if exclude_defocus) using ZernikeCoeff2PTT2
    IrisAO_PTT_res_2          = IrisAO_PTT_res;
    IrisAO_PTT_res_2_residual = IrisAO_PTT - IrisAO_PTT_GT; 
    if exclude_defocus == 1, zex = [1:3, 5]; else, zex = 1:3; end
    for i = zex
        IrisAO_PTT_res_2 = IrisAO_PTT_res_2 + squeeze(zernikeCoeff1(i)/2 * PTT_mode_base(i,:,:));
    end
    IrisAO_PTT_2 = IrisAO_PTT0_2 + IrisAO_PTT_res_2;

    % Regularize locked/weak for nicer display
    IrisAO_PTT_res_2 = CalcLockedSegShift(IrisAO_PTT_res_2, Param.LockedSeg, Param, Convert_hex2jk, Convert_jk2hex);
    IrisAO_PTT_res_2 = CalcLockedSegShift(IrisAO_PTT_res_2, Param.WeakSeg,  Param, Convert_hex2jk, Convert_jk2hex);
    IrisAO_PTT_2     = IrisAO_PTT0_2 + IrisAO_PTT_res_2;

    % Display (round 2)
    wavefront1         = -2 * IrisAO_PTT_res_2(:, 1);
    wavefront2         = -2 * IrisAO_PTT_2(:,   1);
    wavefront_residual = -2 * (IrisAO_PTT_2(:, 1) - IrisAO_PTT_GT(:, 1));
    wavefront_gt       = -2 * IrisAO_PTT_GT(:, 1);

    res_slopes        = IrisAO_PTT_res_2(:,3:-1:2) * 2 / 1000 * SHcam_pitch;
    collective_slopes = IrisAO_PTT_2(:, 3:-1:2)    * 2 / 1000 * SHcam_pitch;
    residual_slopes   = (IrisAO_PTT_2(:,3:-1:2) - IrisAO_PTT_GT(:,3:-1:2)) * 2 / 1000 * SHcam_pitch;
    gt_slopes         = IrisAO_PTT_GT(:,3:-1:2)    * 2 / 1000 * SHcam_pitch;

    wavefront1_d        = weight_display(DM_im, PosRef, wavefront1,        DM_period, Param.UnreachableSeg, res_slopes);
    wavefront2_d        = weight_display(DM_im, PosRef, wavefront2,        DM_period, Param.UnreachableSeg, collective_slopes);
    wavefrontResidual_d = weight_display(DM_im, PosRef, wavefront_residual,DM_period, Param.UnreachableSeg, residual_slopes);
    wavefrontGT_d       = weight_display(DM_im, PosRef, wavefront_gt,      DM_period, Param.UnreachableSeg, gt_slopes);
   
    % Calculate RMS error in wave
    res2_zero = wavefrontResidual_d - min(wavefrontResidual_d(:));
    mask = poly2mask(centerPos(1)+xp, centerPos(2)+yp, size(res2_zero,1), size(res2_zero,2));
    res_vals            = wavefrontResidual_d(mask);
    rms_error2          = sqrt(sum(res_vals.^2)/numel(res_vals)); 
    gt_vals             = wavefrontGT_d(mask);
    gt_rms_error2       = sqrt(sum(gt_vals.^2)/numel(gt_vals));   
    gt_rms_error2_lambda= gt_rms_error2 / 0.515;                   

    figure(6); set(gcf,'Position',[810 450 1100 300]);
    subplot(1,4,1); imagesc(wavefront1_d - min(wavefront1_d(:))); axis image; xticks([]); yticks([]);colorbar;colormap jet;  title('Measured Aberration (post-PTTD)');
    subplot(1,4,2); imagesc(wavefront2_d - min(wavefront2_d(:))); axis image; colorbar;colormap jet;  xticks([]); yticks([]);title('Corrective Wavefront (post-PTTD)');
    subplot(1,4,3); imagesc(wavefrontResidual_d - min(wavefrontResidual_d(:))); axis image; xticks([]); yticks([]); colorbar;colormap jet;  title(['Residual Wavefront: ', num2str(rms_error2)]);

    % Phase wrapping (round 2) for display
    phase_wrap_2 = IrisAO_PTT_2(:, 1);
    while any(phase_wrap_2 >  phase_thres), phase_wrap_2(phase_wrap_2 >  phase_thres) = phase_wrap_2(phase_wrap_2 >  phase_thres) - phase_thres; end
    while any(phase_wrap_2 < -phase_thres), phase_wrap_2(phase_wrap_2 < -phase_thres) = phase_wrap_2(phase_wrap_2 < -phase_thres) + phase_thres; end

    wavefront3_d = weight_display(DM_im, PosRef, -2*(phase_wrap_2 - IrisAO_PTT_GT(:,1)), DM_period, Param.UnreachableSeg, res_slopes); 
    subplot(1,4,4); imagesc(wavefrontGT_d - min(wavefrontGT_d(:))); axis image; colorbar; title('GT Wavefront');
end

%% ------------------ Wrap and zero out outliers; save results --------------
IrisAO_PTT_wrap = IrisAO_PTT;  IrisAO_PTT_wrap(:, 1) = phase_wrap;
if Param.isRmvPTTD
    IrisAO_PTT_2_wrap = IrisAO_PTT_2;  
    IrisAO_PTT_2_wrap(:, 1) = phase_wrap_2;
end

IrisAO_PTT_DM_wrap   = IrisAO_PTT_wrap;      IrisAO_PTT_DM_wrap(Param.Outlier_hex,:)   = 0;
IrisAO_PTT_2_DM_wrap = IrisAO_PTT_2_wrap;    IrisAO_PTT_2_DM_wrap(Param.Outlier_hex,:) = 0;

% Save results locally (portable)
save_name = fullfile(dataDir, 'IrisAO_PTT_2_DM_wrap_result.mat');
save(save_name, "IrisAO_PTT_2_DM_wrap");
save(save_name, 'IrisAO_PTT_2_wrap', 'IrisAO_PTT_2', '-append');
