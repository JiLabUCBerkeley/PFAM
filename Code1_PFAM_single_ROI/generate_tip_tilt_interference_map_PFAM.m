%% PFAM single-ROI processing 
% Purpose:
%   - Preload all (tip, tilt) image stacks 
%   - Build a small ROI mask from weighted centroids in the reference image
%   - Compute per-frame intensity metric inside ROI 
%   - FFT per case, collect normalized/non-normalized peaks at assigned frequencies
%   - Estimate tip/tilt peak positions with a robust threshold + regionprops
%
% Notes:
%   - All optional plots are controlled by DEBUG_PLOTS

close all; clear; clc;

%% -------------------------- Configuration ---------------------------------
DEBUG_PLOTS   = true;         % Toggle diagnostic figures

%% ------------------------- Folder / Path Setup ----------------------------
currentFolder = pwd;
[~, folderName] = fileparts(currentFolder);

% Enforce script location for reproducibility
if ~strcmp(folderName, 'Code1_PFAM_single_ROI')
    error('Run this script from "Code1_PFAM_single_ROI". Current folder: %s', folderName);
end

rootFolder  = fileparts(currentFolder);
funcFolder  = fullfile(rootFolder, 'function');
addpath(funcFolder);

% ---------------------- Choose round folder here ---------------------------
save_folder = fullfile(currentFolder, 'data', 'round3');

% ---------------------- Target segment auto-selection ----------------------
if contains(save_folder, 'round1')
    target_segment = [1 20 23 26 29 32 35 92 95 98 101 104 107 110 113 116 119 122 125];
elseif contains(save_folder, 'round2')
    target_segment = [9 13 17 44 52 60 63 66 73 76 83 86 137 140 151 154 165 168];
elseif contains(save_folder, 'round3')
    target_segment = [11 15 19 40 48 56 68 71 78 81 88 91 130 133 144 147 158 161];
else
    error('Unknown round in save_folder: %s', save_folder);
end
[~, roundName] = fileparts(save_folder);   % extracts "round1", "round2", ...
OUT_MAT_NAME = [roundName, '.mat'];        % e.g. "round2.mat"

% Load indices and pre-allocated sizes
tip_ind_saved   = load(fullfile(save_folder, 'tip_ind_saved.mat')).tip_ind_saved;
tilt_ind_saved  = load(fullfile(save_folder, 'tilt_ind_saved.mat')).tilt_ind_saved;
total_ind_saved = load(fullfile(save_folder, 'total_ind_saved.mat')).total_ind_saved;

n_tip_art  = size(tip_ind_saved,  2);
n_tilt_art = size(tilt_ind_saved, 2);
nCases     = numel(total_ind_saved);
N_REPEAT      = 200;           % Number of frames to keep for each time-series

% Common path prefix for each (tip,tilt) case; each case is suffixed by (tip_idx-1)*n_tip + tilt_idx
case_prefix = fullfile(save_folder, 'Periodp16_ExposureTimep008_nFrames1_length_1000_7Segment_19_power_5_');

%% ---------------------- Preload image stacks (I/O once) -------------------
Image_ref_all = cell(nCases, 1);   % aberrated_Image.tif (reference)
Image_all     = cell(nCases, 1);   % test.tiff (time series)

for ii = 1:nCases
    [tip_idx, tilt_idx] = ind2sub([n_tip_art, n_tip_art], total_ind_saved(ii));
    case_id   = (tip_idx-1)*n_tip_art + tilt_idx;
    case_path = [case_prefix, num2str(case_id)]; 

    Image_ref_all{ii} = double(loadtiff(fullfile(case_path, 'aberrated_Image.tif')));
    Image_all{ii}     = double(loadtiff(fullfile(case_path, 'test.tiff')));
end

%% ----------------------------- ROI builder --------------------------------
% Build ROI from a single reference image (per case), then compute time series.

Intensity_saved = zeros(n_tip_art, n_tilt_art, N_REPEAT);

for ii = 1:nCases
    [tip_idx, tilt_idx] = ind2sub([n_tip_art, n_tip_art], total_ind_saved(ii));
    t0 = tic;

    % Extract stacks for this case
    Iref = Image_ref_all{ii};
    I    = Image_all{ii};

    % Crop borders (2 px all around) and keep first z for ref; all frames for I
    Iref = Iref(2:end-2, 2:end-2, 1);
    I    = I(   2:end-2, 2:end-2, :);
    THR_K = 5;             % Threshold multiplier (mean + K*std) for ROI building

    % Build ROI from weighted centroids of bright blobs
    thr = mean(Iref(:)) + THR_K * std(Iref(:));
    BW  = Iref > thr;
    BW2 = bwareafilt(BW, 5);

    props = regionprops(BW2, Iref, 'WeightedCentroid');
    [H, W] = size(Iref);
    ROI    = false(H, W);
    CROP_SZ       = 2;             % Half-size for square ROI around each centroid (5x5 total)

    for k = 1:numel(props)
        r = max(1, floor(props(k).WeightedCentroid(2)) - CROP_SZ) : ...
            min(H, floor(props(k).WeightedCentroid(2)) + CROP_SZ);
        c = max(1, floor(props(k).WeightedCentroid(1)) - CROP_SZ) : ...
            min(W, floor(props(k).WeightedCentroid(1)) + CROP_SZ);
        ROI(r, c) = true;
    end

    ROI_bg   = ~ROI;
    roi_npix = max(nnz(ROI),    1);
    bg_npix  = max(nnz(ROI_bg), 1);

    % Vectorized per-frame sums
    roi_sum = squeeze(sum(sum(I .* ROI,    1), 2));  % [T x 1]
    bg_sum  = squeeze(sum(sum(I .* ROI_bg, 1), 2));

    % Choose metric (use ROI mean; comment-in background subtraction if desired)
    Metric     = roi_sum / roi_npix;
    % Background = bg_sum  / bg_npix;
    % Intensity  = Metric - Background;
    Intensity  = Metric;

    % Fit to N_REPEAT frames (truncate/pad)
    if numel(Intensity) < N_REPEAT
        Intensity = [Intensity(:); zeros(N_REPEAT - numel(Intensity), 1)];
    else
        Intensity = Intensity(1:N_REPEAT);
    end

    Intensity_saved(tip_idx, tilt_idx, :) = Intensity;

    if DEBUG_PLOTS
        fprintf('(%d/%d) tip=%d, tilt=%d, dt=%.3fs\n', ii, nCases, tip_idx, tilt_idx, toc(t0));       
        figure(200); clf;
        subplot(1,2,1);
        imagesc(Iref); axis image; colormap hot; colorbar;
        title(sprintf('Original Ref (tip=%d, tilt=%d)', tip_idx, tilt_idx));

        subplot(1,2,2);
        imagesc(ROI); axis image; colormap gray; colorbar;
        title('ROI Mask');
        drawnow;

    end
end

save(fullfile(save_folder, 'Intensity_saved.mat'), 'Intensity_saved');

%% ---------------------- FFT analysis / frequency picks --------------------
% Load back 
Intensity_saved = load(fullfile(save_folder, 'Intensity_saved.mat')).Intensity_saved;
Intensity_saved = Intensity_saved(:, :, 3:N_REPEAT);  % discard first two frames for DM segments settlement

% Target segments + frequency list
f = (20:100).';     % candidate frequency bins
f = f(1:length(target_segment));

n_segment = numel(target_segment);

saved_fft_nonnormalized = zeros(n_tip_art, n_tilt_art, n_segment);
saved_fft_normalized    = zeros(n_tip_art, n_tilt_art, n_segment);
save_dc                 = zeros(n_tip_art, n_tilt_art);

for idx = 1:nCases
    [tip_idx, tilt_idx] = ind2sub([n_tip_art, n_tip_art], total_ind_saved(idx));

    y = squeeze(Intensity_saved(tip_idx, tilt_idx, :));
    Y = abs(fft(y));

    % Store normalized and raw magnitudes at target frequencies
    for n = 1:n_segment
        r_tip  = tip_ind_saved(n,  tip_idx);
        r_tilt = tilt_ind_saved(n, tilt_idx);
        bin    = f(n) + 1;               % MATLAB 1-based index for FFT bins

        saved_fft_nonnormalized(r_tip, r_tilt, n) = Y(bin);
        saved_fft_normalized(   r_tip, r_tilt, n) = Y(bin) / max(Y(1), eps);
    end
end

% Save maps if desired later:
save(fullfile(save_folder, 'tip_tilt_map.mat'), 'saved_fft_nonnormalized');
%% ---------------------- Peak localization over segments -------------------
% Locate peaks on each segment’s non-normalized FFT map and convert to angle.
pos_ref   = [ceil(n_tip_art/2), ceil(n_tilt_art/2)];
pix_angle = 0.15 / 0.6306;  % (example) pixel-to-angle conversion (rad/pixel or similar)

tip_tilt_saved = zeros(n_segment, 2);
shift_saved    = zeros(n_segment, 2);

for ii = 1:n_segment
    M = saved_fft_nonnormalized(:, :, ii);

    % Normalize to [0, 1] for robust thresholding
    z = (M - min(M(:))) / max(max(M(:)) - min(M(:)), eps);

    thr = mean(z(:)) + 1 * std(z(:));  % mild threshold for peak regions
    BW  = z > thr;

    % Keep labeled blobs within a reasonable size and prioritize high-intensity ones
    BW  = bwareafilt(BW, [1, 80], 4);                % size constraint (px)
    BW  = bwpropfilt(BW, z, 'MaxIntensity', 1, "largest", 4);

    props = regionprops(BW, z, 'WeightedCentroid');

    if isempty(props)
        % No reliable blob found; mark as NaN
        shift_saved(ii, :) = [NaN, NaN];
    else
        % Use the strongest blob (props(1) after bwpropfilt)
        c = props(1).WeightedCentroid;          % [x, y] in MATLAB coords
        shift_saved(ii, :) = [c(2), c(1)] - pos_ref;  % [rowShift, colShift]
    end

    tip_tilt_saved(ii, :) = shift_saved(ii, :) * pix_angle;

    if DEBUG_PLOTS
        figure(10); subplot(4,5,ii);
        imagesc(z); axis image; colormap jet; colorbar; title(sprintf('Seg %d', ii)); hold on;
        if all(isfinite(shift_saved(ii, :)))
            plot(shift_saved(ii,2)+pos_ref(2), shift_saved(ii,1)+pos_ref(1), 'kx', 'MarkerSize', 12, 'LineWidth', 2);
        end
        hold off
    end
end

%% ------------------------------- Save -------------------------------------
save(fullfile(save_folder, OUT_MAT_NAME), 'tip_tilt_saved');

if DEBUG_PLOTS
    fprintf('Saved tip/tilt to: %s\n', fullfile(save_folder, OUT_MAT_NAME));
end
