%% PFAM SIFT multiple-ROIs
% Purpose
%   1) Build ROIs from a reference image (two sub-windows shown below)
%   2) For each (tip, tilt) case, compute per-frame FT power spectrum metric
%      by summing six sideband crops in the frequency domain
%   3) FFT the time series per case; assemble per-segment maps
%   4) Localize peaks to produce per-segment (tip, tilt) estimates
%
% Notes
%   - Plots are gated by DEBUG_PLOTS

close all; clear; clc;

%% ------------------------------ Configuration ------------------------------
DEBUG_PLOTS   = true;         % Toggle diagnostic figures
ROUND_NAME    = 'round2';     % 'round1' | 'round2'
N_REPEAT      = 200;          % Frames to keep from test.tiff
CROP_FT       = 3;            % Half-width of FT crop window (sideband sum)
PERIOD_BASE   = 5.95;         % Lattice period for multiple-ROI workflow

% Optical geometry for sideband radius
PIX_PITCH_CAM = 6.5;                            % μm
MAG_SLM_CAM   = (125/150) * (300/175) * (75/85);

%% ------------------------------ Paths & Guards -----------------------------
here = pwd;
[~, folderName] = fileparts(here);
if ~strcmp(folderName, 'Code4_PFAM_SIFT_multiple_ROIs')
    error('Run this script from "Code4_PFAM_SIFT_multiple_ROIs". Current folder: %s', here);
end

dataFolder = fullfile(here, 'data');
saveFolder = fullfile(dataFolder, ROUND_NAME);
if ~exist(saveFolder, 'dir')
    error('Missing data folder: %s', saveFolder);
end

% Index files
tip_ind_path   = fullfile(saveFolder, 'tip_ind_saved.mat');
tilt_ind_path  = fullfile(saveFolder, 'tilt_ind_saved.mat');
total_ind_path = fullfile(saveFolder, 'total_ind_saved.mat');
assert(all([exist(tip_ind_path,'file'), exist(tilt_ind_path,'file'), exist(total_ind_path,'file')] == 2), ...
    'Missing one or more index files in %s', saveFolder);

tip_ind_saved   = load(tip_ind_path).tip_ind_saved;
tilt_ind_saved  = load(tilt_ind_path).tilt_ind_saved;
total_ind_saved = load(total_ind_path).total_ind_saved;

n_tip_art  = size(tip_ind_saved,  2);
n_tilt_art = size(tilt_ind_saved, 2);
nCases     = numel(total_ind_saved);

% Case path prefix (e.g., ..._5_17 for tip=5, tilt=17)
case_prefix = fullfile(saveFolder, 'Periodp16_ExposureTimep008_nFrames1_length_1000_7Segment_19_power_5_');

%% ------------------------------ Load Template & Sub-ROIs --------------------
% Use case #1 reference for ROI planning
case1_path   = [case_prefix, '1'];
refPath_case = fullfile(case1_path, 'reference_Image.tif');
if ~exist(refPath_case, 'file')
    refPath_case = fullfile(case1_path, 'aberrated_Image.tif');  % fallback
end
assert(exist(refPath_case,'file')==2, 'Missing reference image in %s', case1_path);

Iref_full = double(loadtiff(refPath_case));
Iref_full = Iref_full(4:end-5, 4:end-5);    % Border crop 

% Choose two sub-ROIs (edit here if needed)
sub_num  = 2;
subrange = zeros(sub_num, 4);     % [r1 r2 c1 c2]
crop1    = 32;
crop2    = 32;

% (Example) You can add different subranges for own data
if strcmpi(ROUND_NAME, 'round1')
    subrange(1,:) = [106 - crop1 + 1, 106 + crop1, 189 - crop1 + 1, 189 + crop1];
    subrange(2,:) = [128 - crop2 + 1, 128 + crop2, 1379 - crop2 + 1, 1379 + crop2];
else
    subrange(1,:) = [106-10 - crop1 + 1, 106-10 + crop1, 189-10 - crop1 + 1, 189-10 + crop1];
    subrange(2,:) = [128-10 - crop2 + 1, 128-10 + crop2, 1379-10 - crop2 + 1, 1379-10 + crop2];
end

if DEBUG_PLOTS
    figure(13310);
    imagesc(Iref_full); axis image; colormap hot; title('Reference (cropped) with sub-ROIs');
    hold on;
    for i = 1:sub_num
        r1 = subrange(i,1); r2 = subrange(i,2);
        c1 = subrange(i,3); c2 = subrange(i,4);
        rectangle('Position', [c1, r1, c2-c1+1, r2-r1+1], 'EdgeColor','w');
    end
    hold off;
end

%% ------------------------------ Inspect each sub-ROI for example image ------------------------
% Build binary masks + display FT crops for sanity
for i = 1:sub_num
    r1 = subrange(i,1); r2 = subrange(i,2);
    c1 = subrange(i,3); c2 = subrange(i,4);

    Iroi = medfilt2(Iref_full(r1:r2, c1:c2), [3,3]);

    if DEBUG_PLOTS
        figure(13311); subplot(sub_num,1,i); imagesc(Iroi); axis image; colormap hot; title(sprintf('Sub-ROI %d', i));
    end

    % SI Frequency locations (6 sidebands)
    H = size(Iroi,1); W = size(Iroi,2);
    ft_img = fftshift(fft2(fftshift(Iroi)));
    A2     = abs(ft_img).^2;                                     % (no square — visualization only)

    period_slm = PERIOD_BASE * 8.2;
    rho        = (period_slm * MAG_SLM_CAM)^(-1) / ( (H * PIX_PITCH_CAM)^(-1) ) * ones(4,1);
    ang        = zeros(4,1);
    for j = 0:2
        ang(j+1) = atan2(12,-1) - j*2*pi/3;
    end
    ang(4) = ang(1);

    xoff = round(rho .* cos(ang));
    yoff = round(rho .* sin(ang));
    cx   = floor(H/2) + 1; cy = floor(W/2) + 1;

    cens = [ ...
        cx - yoff(1), cy + xoff(1);
        cx + yoff(1), cy - xoff(1);
        cx - yoff(2), cy + xoff(2);
        cx + yoff(2), cy - xoff(2);
        cx - yoff(3), cy + xoff(3);
        cx + yoff(3), cy - xoff(3) ];

    if DEBUG_PLOTS
        figure(122); clf; imagesc(A2); axis image; colormap jet; colorbar; title(sprintf('FT (log) sub-ROI %d', i));
        set(gca,'ColorScale','log'); hold on; plot(cens(:,2), cens(:,1), 'b*'); hold off;

        % Show 6 cropped windows
        figure(123); clf;
        for ci = 1:6
            rr = max(cens(ci,1)-CROP_FT,1):min(cens(ci,1)+CROP_FT,H);
            cc = max(cens(ci,2)-CROP_FT,1):min(cens(ci,2)+CROP_FT,W);
            subplot(3,2,ci); imagesc(A2(rr,cc)); axis image; colormap jet; colorbar;
            set(gca,'ColorScale','log'); title(sprintf('Crop %d', ci));
        end
    end
end

%% ------------------------------ Metric per Case / ROI (FAST) ---------------
% Vectorized over frames:
%   - One fft2 per stack (applied page-wise)
%   - One masked sum per ROI (six sidebands union mask)
Intensity_saved = zeros(sub_num, n_tip_art, n_tilt_art, N_REPEAT, 'double');

% --- Precompute sideband masks per sub-ROI (size depends only on ROI shape)
roiMasks = cell(sub_num,1);
for s = 1:sub_num
    r1 = subrange(s,1); r2 = subrange(s,2);
    c1 = subrange(s,3); c2 = subrange(s,4);
    Hs = r2 - r1 + 1;  Ws = c2 - c1 + 1;
    cx = floor(Hs/2) + 1;  cy = floor(Ws/2) + 1;

    % radial offsets (six points across three axes ±)
    period_slm = PERIOD_BASE * 8.2;
    rho  = (period_slm * MAG_SLM_CAM)^(-1) / ((Hs * PIX_PITCH_CAM)^(-1)) * ones(4,1);
    ang  = zeros(4,1);
    for j = 0:2, ang(j+1) = atan2(12,-1) - j*2*pi/3; end
    ang(4) = ang(1);
    xoff = round(rho .* cos(ang));
    yoff = round(rho .* sin(ang));

    centers = [ ...
        cx - yoff(1), cy + xoff(1);
        cx + yoff(1), cy - xoff(1);
        cx - yoff(2), cy + xoff(2);
        cx + yoff(2), cy - xoff(2);
        cx - yoff(3), cy + xoff(3);
        cx + yoff(3), cy - xoff(3) ];

    % Build a single binary mask that is the union of six CROP_FT squares
    M = false(Hs, Ws);
    for ci = 1:6
        rr = max(centers(ci,1)-CROP_FT,1):min(centers(ci,1)+CROP_FT,Hs);
        cc = max(centers(ci,2)-CROP_FT,1):min(centers(ci,2)+CROP_FT,Ws);
        M(rr,cc) = true;
    end
    roiMasks{s} = M;
end

% --- Loop (or parfor) over cases; inside is fully vectorized over frames
% parfor idx = 1:nCases   % <- enable if you have Parallel Toolbox
for idx = 1:nCases
    [tip_idx, tilt_idx] = ind2sub([n_tip_art, n_tilt_art], total_ind_saved(idx));
    case_id   = (tip_idx-1)*n_tip_art + tilt_idx;
    case_path = [case_prefix, num2str(case_id)];

    % Load stacks once
    Iref = double(loadtiff(fullfile(case_path,'aberrated_Image.tif')));
    I    = double(loadtiff(fullfile(case_path,'test.tiff')));

    % Border crop
    Iref = Iref(4:end-5, 4:end-5, 1);
    I    = I(   4:end-5, 4:end-5, :);

    for s = 1:sub_num
        r1 = subrange(s,1); r2 = subrange(s,2);
        c1 = subrange(s,3); c2 = subrange(s,4);

        % Extract ROI stack (Hs x Ws x T)
        Iroi = I(r1:r2, c1:c2, :);

        % Page-wise 2D FFT on all frames at once, with pre/post fftshift on dims 1&2
        % (Use two calls to fftshift to specify dims explicitly)
        Iroi_c = fftshift(fftshift(Iroi,1),2);
        F      = fft2(Iroi_c);                             % FFT per page
        A2     = abs(fftshift(fftshift(F,1),2)).^2;        % centered power spectrum, per page

        % Masked sum over (row,col), keep time dimension
        M = roiMasks{s};
        % Expand mask to 3D without copying (implicit expansion), then sum over dims 1&2
        S = sum(A2 .* M, [1 2]);                           % 1x1xT
        Metric = squeeze(S);                                % Tx1

        % Truncate/pad to N_REPEAT
        if numel(Metric) < N_REPEAT
            Metric(N_REPEAT) = 0;
        else
            Metric = Metric(1:N_REPEAT);
        end

        Intensity_saved(s, tip_idx, tilt_idx, :) = Metric;
    end

    if DEBUG_PLOTS
        T = size(I,3);
        fprintf('Case (%3d/%3d): tip=%2d, tilt=%2d, frames=%d\n', idx, nCases, tip_idx, tilt_idx, T);
    end
end

save(fullfile(saveFolder, 'Intensity_saved_multiROI.mat'), 'Intensity_saved');

%% ------------------------------ FFT → maps ----------------------------------
% Use frames 3:N to allow settling
Intensity_saved = Intensity_saved(:,:,:,3:N_REPEAT);

% Target segments per round
switch lower(ROUND_NAME)
    case 'round1'
        target_segment = [1 20 23 26 29 32 35];
    case 'round2'
        target_segment = [92 95 98 101 104 107 110 113 116 119 122 125];
end

f = (20:100)';                              % candidate frequency bins
f = f(1:length(target_segment));

n_segment = numel(target_segment);
saved_fft_nonnormalized = zeros(sub_num, n_tip_art, n_tilt_art, n_segment);
saved_fft_normalized    = zeros(sub_num, n_tip_art, n_tilt_art, n_segment);
save_dc                 = zeros(n_tip_art, n_tilt_art);

for s = 1:sub_num
    Is = squeeze(Intensity_saved(s, :, :, :));   % [tip x tilt x T]

    for idx = 1:nCases
        [tip_idx, tilt_idx] = ind2sub([n_tip_art, n_tilt_art], total_ind_saved(idx));

        y = squeeze(Is(tip_idx, tilt_idx, :));
        Y = abs(fft(y));

        for n = 1:n_segment
            r_tip  = tip_ind_saved(n,  tip_idx);
            r_tilt = tilt_ind_saved(n, tilt_idx);
            bin    = f(n) + 1;

            saved_fft_nonnormalized(s, r_tip, r_tilt, n) = Y(bin);
            saved_fft_normalized(   s, r_tip, r_tilt, n) = Y(bin) / max(Y(1), eps);
        end

        save_dc(tip_idx, tilt_idx) = Y(1);
    end
end

save(fullfile(saveFolder, 'tip_tilt_map_multiROI.mat'), 'saved_fft_nonnormalized');

if DEBUG_PLOTS
    roi_to_show = 1;
    figure(12311);
    for n = 1:n_segment
        subplot(4,5,n);
        imagesc(medfilt2(squeeze(saved_fft_normalized(roi_to_show,:,:,n))));
        axis image; colormap gray; title(sprintf('ROI%d Seg %d', roi_to_show, n)); colorbar;
    end
    figure(12312);
    for n = 1:n_segment
        subplot(4,5,n);
        imagesc(medfilt2(squeeze(saved_fft_nonnormalized(roi_to_show,:,:,n))));
        axis image; colormap gray; title(sprintf('ROI%d Seg %d', roi_to_show, n)); colorbar;
    end
    figure(12313); imagesc(save_dc); axis image; colormap gray; colorbar; title('DC map');
end

%% ------------------------------ Peak localization ---------------------------
% For each ROI, localize the tip/tilt peak on the non-normalized FFT map
PIX_TO_ANGLE  = 0.15 / 0.6306;
pos_ref = [ceil(n_tip_art/2), ceil(n_tilt_art/2)];
[Xg, Yg] = meshgrid(-4:4, -4:4); XY(:,:,1) = Xg; XY(:,:,2) = Yg;
gaussfn  = @(B,XY) B(5) * exp(-((XY(:,:,1)-B(2))/B(3)).^2 - ((XY(:,:,2)-B(1))/B(4)).^2);

for s = 1:sub_num
% for s = 1

    tip_tilt_saved = zeros(n_segment, 2);

    resnorm_save = zeros(1, n_segment);
    for n = 1:n_segment
    % for n = 1

        M = squeeze(saved_fft_nonnormalized(s,:,:,n));
        z = (M - min(M(:))) / max(max(M(:)) - min(M(:)), eps);
        thr  = mean(z(:)) + 1*std(z(:));
        BW   = z > thr;
        BW   = bwareafilt(BW, [4, 80], 4);
        BW   = bwpropfilt(BW, z, 'MaxIntensity', 1, "largest", 4);

        props = regionprops(BW, z, 'WeightedCentroid');

        if isempty(props)
            % Fallback to Gaussian fit if no blob found
            B0 = [0, 0, 0.5, 0.5, 1];
            lb = [-4, -4, 0,   0,   0];
            ub = [ 4,  4, 4,   4,  10];
            try
                [B, resn] = lsqcurvefit(gaussfn, B0, XY, z, lb, ub);
            catch
                B = [0 0 0 0 0]; resn = inf;
            end
            resnorm_save(n) = resn;
            if resn > 2
                tip_tilt_saved(n,:) = [NaN, NaN];
            else
                shift = B(1:2); % [dy, dx]
                tip_tilt_saved(n,:) = shift * PIX_TO_ANGLE;
            end
        else
            c      = props(1).WeightedCentroid;         % [x, y]
            shift  = [c(2), c(1)] - pos_ref;            % [dy, dx]
            tip_tilt_saved(n,:) = shift * PIX_TO_ANGLE;
        end

        if DEBUG_PLOTS
            figure(12314);
            subplot(4,5,n); imagesc(M); axis image; colormap jet; colorbar;
            title(sprintf('ROI%d Seg %d', s, n)); hold on;
            if all(isfinite(tip_tilt_saved(n,:)))
                plot(shift(2)+pos_ref(2), shift(1)+pos_ref(1), 'kx', 'MarkerSize', 10, 'LineWidth', 2);
            end
            hold off

            figure(12315);
            subplot(4,5,n); imagesc(BW); axis image; colormap gray; colorbar;
            title(sprintf('Mask ROI%d Seg %d', s, n));
        end
    end

    if DEBUG_PLOTS
        figure(123118); plot(resnorm_save, '-o'); grid on;
        title(sprintf('ROI %d Gaussian fit residuals', s));
        xlabel('Segment'); ylabel('resnorm');
    end

    % Save per-ROI result
    outName = sprintf('%s_centroid_multiROI_%02d.mat', ROUND_NAME, s);
    save(fullfile(saveFolder, outName), 'tip_tilt_saved');
    fprintf('Saved: %s\n', fullfile(saveFolder, outName));
end
