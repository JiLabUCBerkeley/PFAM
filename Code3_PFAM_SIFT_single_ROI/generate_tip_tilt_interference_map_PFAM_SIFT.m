%% PFAM SIFT single-ROI processing
% Purpose:
%   - Build a bright-blob ROI from a reference image
%   - Construct 6 FT sideband windows from lattice geometry
%   - For each (tip, tilt) case, compute per-frame Power spectrum metric
%   - FFT the temporal metric and assemble tip/tilt maps
%   - Localize peak per target segment to produce (tip, tilt) in angle units
%
% Notes:
%   - All plots are gated by DEBUG_PLOTS

close all; clear; clc;

%% ------------------------------ Configuration ------------------------------
DEBUG_PLOTS = true;                 % Toggle figures
ROUND_NAME  = 'round1';             % 'round1' | 'round2' | 'round3'
N_REPEAT    = 200;                  % Max frames per time series (test.tiff)
CROP_FT     = 3;                    % Half width of crop window on FT peaks
THR_K       = 1.0;                  % z-threshold multiplier for blob (mean + K*std)


% (Lattice) period parameters for sideband locations (leave numeric as-is)
period                = 3.97;       % pixel period at slm
period_slm            = period * 8.2; % slm pixel size: 8.2um
mag_Slm_Camera        = (125/150) * (300/175) * (75/85);

%% ------------------------------ Path & Guards ------------------------------
currentFolder = pwd;
[~, folderName] = fileparts(currentFolder);
if ~strcmp(folderName, 'Code3_PFAM_SIFT_single_ROI')
    error('Run this script from "Code3_PFAM_SIFT_single_ROI". Current folder: %s', folderName);
end

dataRoot   = fullfile(currentFolder, 'data');
save_root  = fullfile(dataRoot, ROUND_NAME);
if ~exist(save_root, 'dir'), error('Missing data folder: %s', save_root); end

% Files required
tipFile   = fullfile(save_root, 'tip_ind_saved.mat');
tiltFile  = fullfile(save_root, 'tilt_ind_saved.mat');
totalFile = fullfile(save_root, 'total_ind_saved.mat');
assert(exist(tipFile,'file')==2 && exist(tiltFile,'file')==2 && exist(totalFile,'file')==2, ...
    'Missing tip/tilt/total index files in %s', save_root);

% Load index maps
tip_ind_saved   = load(tipFile).tip_ind_saved;
tilt_ind_saved  = load(tiltFile).tilt_ind_saved;
total_ind_saved = load(totalFile).total_ind_saved;

n_tip_art  = size(tip_ind_saved,  2);
n_tilt_art = size(tilt_ind_saved, 2);
nCases     = numel(total_ind_saved);

% Case path prefix (e.g. ..._5_17 for tip=5, tilt=17)
case_prefix = fullfile(save_root, 'Periodp16_ExposureTimep008_nFrames1_length_1000_7Segment_19_power_5_');

%% ------------------------------ Template & ROI -----------------------------
% Use case #1 for template/reference
case1_path   = [case_prefix, '1'];
refPath_case = fullfile(case1_path, 'reference_Image.tif');
if ~exist(refPath_case, 'file')
    % Some datasets use 'aberrated_Image.tif' for the static reference
    refPath_case = fullfile(case1_path, 'aberrated_Image.tif');
    assert(exist(refPath_case,'file')==2, 'Missing both reference_Image.tif and aberrated_Image.tif in %s', case1_path);
end

Iref_full = double(loadtiff(refPath_case));
Iref      = Iref_full(2:end-2, 2:end-2);         % border crop


% Precompute FT geometry for 6 sidebands
ft_img    = fftshift(fft2(fftshift(Iref)));
abs_ft    = abs(ft_img).^2; % Calculate Power spectrum
H         = size(Iref,1);
W         = size(Iref,2);
cx        = floor(H/2) + 1;
cy        = floor(W/2) + 1;

rho   = ( (period_slm*mag_Slm_Camera)^-1 ) / ( (H*6.5)^-1 ) * ones(4,1);  % 6.5um: camera pixel pitch
angles = zeros(4,1);
for j = 0:2
    angles(j+1) = atan2(12, -1) - j*2*pi/3;
end
angles(4) = angles(1);

% Convert polar offsets to pixel offsets
x_off = round(rho .* cos(angles));
y_off = round(rho .* sin(angles));

% Sideband centers (6 total, symmetric pairs around DC)
cens = [ ...
   cx - y_off(1), cy + x_off(1);  % 1
   cx + y_off(1), cy - x_off(1);  % 2
   cx - y_off(2), cy + x_off(2);  % 3
   cx + y_off(2), cy - x_off(2);  % 4
   cx - y_off(3), cy + x_off(3);  % 5
   cx + y_off(3), cy - x_off(3)]; % 6
if DEBUG_PLOTS
    figure(12); imagesc(abs_ft); axis image; colormap jet; colorbar; title('Reference FT (log)');
    set(gca,'ColorScale','log'); hold on;
    plot(cens(:,2), cens(:,1), 'b*'); hold off;
end

%% ------------------------------ Per-case Metric ----------------------------
% FT-energy metric per frame (sum over 6 cropped sidebands)
Intensity_saved = zeros(n_tip_art, n_tilt_art, N_REPEAT);

for n = 1:nCases
    [tip_idx, tilt_idx] = ind2sub([n_tip_art, n_tilt_art], total_ind_saved(n));
    case_id   = (tip_idx-1)*n_tip_art + tilt_idx;
    case_path = [case_prefix, num2str(case_id)];

    % Time series (test.tiff). Fall back to existing names as needed.
    testPath = fullfile(case_path, 'test.tiff');
    assert(exist(testPath,'file')==2, 'Missing test.tiff in %s', case_path);

    % Load stacks and crop borders
    I = double(loadtiff(testPath));
    I = I(2:end-2, 2:end-2, :);

    T = size(I,3);
    Metric = zeros(min(T, N_REPEAT), 1);

    for k = 1:min(T, N_REPEAT)
        Ik = I(:,:,k);
        ftk = fftshift(fft2(fftshift(Ik)));
        A2  = abs(ftk).^2;

        % Sum energy in 6 crops
        S = 0;
        for ci = 1:size(cens,1)
            rr = max(cens(ci,1)-CROP_FT,1) : min(cens(ci,1)+CROP_FT,H);
            cc = max(cens(ci,2)-CROP_FT,1) : min(cens(ci,2)+CROP_FT,W);
            S  = S + sum(A2(rr,cc), 'all');
        end

        Metric(k) = S;
    end

    % pad/trim to N_REPEAT
    if numel(Metric) < N_REPEAT
        Metric = [Metric; zeros(N_REPEAT-numel(Metric),1)];
    end

    Intensity_saved(tip_idx, tilt_idx, :) = Metric;

    if DEBUG_PLOTS
        fprintf('(%d/%d) tip=%d tilt=%d: frames=%d\n', n, nCases, tip_idx, tilt_idx, T);
    end
end

save(fullfile(save_root, 'Intensity_saved_FT.mat'), 'Intensity_saved');

%% ------------------------------ FFT & Maps ---------------------------------
% Use frames 3:N to allow settling
Intensity_saved = Intensity_saved(:,:,3:N_REPEAT);
% Select target segments per round (match your convention)
switch lower(ROUND_NAME)
    case 'round1'
        target_segment = [1 20 23 26 29 32 35];
    case 'round2'
        target_segment = [92 95 98 101 104 107 110 113 116 119 122 125];
end

% Frequency picks (bins)
f = (20:100)';                           % candidates
f = f(1:length(target_segment));
n_segment = numel(target_segment);

saved_fft_nonnormalized = zeros(n_tip_art, n_tilt_art, n_segment);
saved_fft_normalized    = zeros(n_tip_art, n_tilt_art, n_segment);
save_dc                 = zeros(n_tip_art, n_tilt_art);

for n = 1:nCases
    [tip_idx, tilt_idx] = ind2sub([n_tip_art, n_tilt_art], total_ind_saved(n));

    y = squeeze(Intensity_saved(tip_idx, tilt_idx, :));
    Y = abs(fft(y));

    for s = 1:n_segment
        r_tip  = tip_ind_saved(s,  tip_idx);
        r_tilt = tilt_ind_saved(s, tilt_idx);
        bin    = f(s) + 1;

        saved_fft_nonnormalized(r_tip, r_tilt, s) = Y(bin);
        saved_fft_normalized(   r_tip, r_tilt, s) = Y(bin) / max(Y(1), eps);
    end

    save_dc(tip_idx, tilt_idx) = Y(1);
end

save(fullfile(save_root, 'tip_tilt_map.mat'), 'saved_fft_nonnormalized');

if DEBUG_PLOTS
    figure(12311);
    for s = 1:n_segment
        subplot(4,5,s); imagesc(medfilt2(saved_fft_normalized(:,:,s)));
        axis image; colormap gray; title(sprintf('Seg %d', s)); colorbar;
    end
    figure(12312);
    for s = 1:n_segment
        subplot(4,5,s); imagesc(medfilt2(saved_fft_nonnormalized(:,:,s)));
        axis image; colormap gray; title(sprintf('Seg %d', s)); colorbar;
    end
end

%% ------------------------------ Peak Localization --------------------------
pos_ref         = [ceil(n_tip_art/2), ceil(n_tilt_art/2)];
tip_tilt_saved  = zeros(n_segment, 2);
shift_saved     = zeros(n_segment, 2);
% Pixel->angle scaling used later when converting shifts (if needed)
PIX_TO_ANGLE = 0.15 / 0.6306;       % example scaling

for s = 1:n_segment
    M = saved_fft_nonnormalized(:,:,s);
    z = (M - min(M(:))) / max(max(M(:)) - min(M(:)), eps);

    thr = mean(z(:)) + 1 * std(z(:));
    BW  = z > thr;
    BW  = bwareafilt(BW, [3, 80], 4);
    BW  = bwpropfilt(BW, z, 'MaxIntensity', 1, 'largest', 4);

    props = regionprops(BW, z, 'WeightedCentroid');

    if isempty(props)
        shift_saved(s,:) = [NaN, NaN];
    else
        c = props(1).WeightedCentroid;            % [x, y]
        shift_saved(s,:) = [c(2), c(1)] - pos_ref; % [dy, dx]
    end

    tip_tilt_saved(s,:) = shift_saved(s,:) * PIX_TO_ANGLE;

    if DEBUG_PLOTS
        figure(10); subplot(4,5,s);
        imagesc(z); axis image; colormap jet; colorbar; title(sprintf('Seg %d', s)); hold on;
        if all(isfinite(shift_saved(s,:)))
            plot(shift_saved(s,2)+pos_ref(2), shift_saved(s,1)+pos_ref(1), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
        end
        hold off
    end
end

%% ------------------------------ Save Output --------------------------------
outName = sprintf('%s_SIFT_ft.mat', ROUND_NAME);
save(fullfile(save_root, outName), 'tip_tilt_saved');
fprintf('Saved: %s\n', fullfile(save_root, outName));
