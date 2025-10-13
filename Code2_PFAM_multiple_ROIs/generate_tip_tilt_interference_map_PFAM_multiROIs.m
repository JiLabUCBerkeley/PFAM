%%
close all; clear; clc;

%% ------------------------------ Configuration --------------------------------
ROUND_NAME            = 'round2';     % choose data round folder round1 / round2
LOAD_PREV_TEMPLATE    = true;         % false for round 1 / true for round 2
DEBUG_PLOT            = false;

%% ------------------------------ Path Setup & Sanity Checks -------------------
currentFolder = pwd;
[~, folderName] = fileparts(currentFolder);
if ~strcmp(folderName, 'Code2_PFAM_multiple_ROIs')
    error('Run this script from "Code2_PFAM_multiple_ROIs". Current folder: %s', folderName);
end

data_Folder       = fullfile(currentFolder, 'data');
save_folder      = fullfile(data_Folder, ROUND_NAME);
templatePrevPath = fullfile(data_Folder, 'template_Image.tif');      % previous template (optional)
beadPosPath      = fullfile(data_Folder, 'bead_position.mat');        % centroids

% Basic existence checks for critical inputs
if ~exist(beadPosPath, 'file')
    error('Missing centroid file: %s', beadPosPath);
end
if ~exist(save_folder, 'dir')
    error('Missing save folder: %s', save_folder);
end

% ---------------------- Target segment auto-selection ----------------------
if contains(save_folder, 'round1')
    target_segment = [1 20 23 26 29 32 35];
elseif contains(save_folder, 'round2')
    target_segment = [92 95 98 101 104 107 110 113 116 119 122 125];
else
    error('Unknown round in save_folder: %s', save_folder);
end

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

%% ------------------------------ Build Template (Sum 81 stacks) --------------
template_Image = 0;
for kk = 1:nCases
    caseFolder = sprintf('Periodp16_ExposureTimep008_nFrames1_length_1000_7Segment_19_power_5_%d', kk);
    tiffPath   = fullfile(save_folder, caseFolder, 'aberrated_Image.tif');
    if ~exist(tiffPath, 'file')
        error('Missing image stack: %s', tiffPath);
    end

    I = double(loadtiff(tiffPath));
    % Normalize to [0,1]
    I = I ./ max(I(:));
    I = I - min(I(:));
    template_Image = template_Image + I;
end

% Re-normalize the summed template to [0,1]
template_Image = template_Image ./ max(template_Image(:));
template_Image = template_Image - min(template_Image(:));
% Suppress Hot Pixels 
MANUAL_HOTPIX_THR = 0.99;
template_median = medfilt2(template_Image, [5 5]);
hot_mask        = template_Image > MANUAL_HOTPIX_THR;
template_Image(hot_mask) = template_median(hot_mask);

% Final normalize
template_Image = template_Image ./ max(template_Image(:));
template_Image = template_Image - min(template_Image(:));

% Save template Image
if ~LOAD_PREV_TEMPLATE 
    imwrite(template_Image, fullfile(data_Folder, 'template_Image.tif'));
end
%% ------------------------------ Load Bead position -------------------------------
S = load(beadPosPath);
if ~isfield(S, 'centroid_save')
    error('Variable "centroid_save" not found in %s', beadPosPath);
end
centroid_save = S.centroid_save;     % Nx2, [x, y]

%% ------------------------------ Estimate Global Shift (Optional) ------------
if LOAD_PREV_TEMPLATE && exist(templatePrevPath, 'file')
    template_prev = double(imread(templatePrevPath));
    template_prev = template_prev ./ max(template_prev(:));
    template_prev = template_prev - min(template_prev(:));

    % Cross-correlation (template_prev is the "template", template_Image is the "image")
    c = normxcorr2(template_prev, template_Image);
    [peak_y, peak_x] = find(c == max(c(:)));
    shift_y = peak_y - size(template_prev, 1);
    shift_x = peak_x - size(template_prev, 2);

    fprintf('Global translation (previous -> current): Δx = %d px, Δy = %d px\n', shift_x, shift_y);

    % Apply shift to centroids
    centroid_shifted = centroid_save + [shift_x, shift_y];
else
    centroid_shifted = centroid_save;
end

%% ------------------------------ Visualize Centroids -------------------------
I_gray = template_Image;                          % shorthand
figure; imagesc(I_gray); axis image; colormap gray; caxis([0, 1]);
title('Template with Centroids'); hold on;
for k = 1:size(centroid_shifted, 1)
    plot(centroid_shifted(k, 1), centroid_shifted(k, 2), 'b*', 'MarkerSize', 8);
end
hold off;

%% ------------------------------ Build Patches -------------------------------
% NOTE: function - divide_image_patches at the end of this code
%       It expects centroids in [row, col] or [y, x]; your prior usage fliplr()
%       indicates it wants [row, col], hence the flip.

PATCH_PAD             = 7;            % half-size padding used by divide_image_patches
PATCH_CLUSTER_DIST    = 25;           % clustering distance (px) for grouping centroids
patches = divide_image_patches(I_gray, fliplr(centroid_shifted), PATCH_CLUSTER_DIST, PATCH_PAD);

% Quick visualization of bounding boxes
figure; imagesc(I_gray); axis image; colormap gray; caxis([0, 0.5]); hold on;
for k = 1:numel(patches)
    bb = patches{k}.bbox;      % [r_min, r_max, c_min, c_max]
    x  = bb(3);                % c_min
    y  = bb(1);                % r_min
    w  = bb(4) - bb(3) + 1;    % width
    h  = bb(2) - bb(1) + 1;    % height
    rectangle('Position', [x, y, w, h], 'EdgeColor', 'r', 'LineWidth', 1.5);
    text(x, y - 5, sprintf('%d', k), 'Color', 'y', 'FontSize', 10, 'FontWeight', 'bold');
end
title('Patches (Bounding Boxes)'); hold off;
if ~LOAD_PREV_TEMPLATE
    % Save Patches 
    patchesSavePath  = fullfile(data_Folder, 'patches.mat');
    save(patchesSavePath, 'patches');
    fprintf('Saved patches to: %s\n', patchesSavePath);
end
%%
num_beads = size(centroid_shifted, 1);      % number of beads (before patch-averaging)

% Storage (per-bead time series)
Intensity_saved = zeros(num_beads, n_tip_art, n_tilt_art, N_REPEAT);

% Common case folder prefix
case_prefix = fullfile(save_folder, 'Periodp16_ExposureTimep008_nFrames1_length_1000_7Segment_19_power_5_');

% --------------------------- Preload stacks to avoid repeated I/O ---------------------
Image_all = cell(nCases, 1);
for k = 1:nCases
    [ti, to] = ind2sub([n_tip_art, n_tip_art], total_ind_saved(k));
    caseFolder = sprintf('%s%d', case_prefix, (ti - 1) * n_tip_art + to);
    stackPath  = fullfile(caseFolder, 'test.tiff');               % test stack
    if ~exist(stackPath, 'file')
        error('Missing image stack: %s', stackPath);
    end
    Image_all{k} = double(loadtiff(stackPath));
end

% --------------------------- Compute per-bead time series ----------------------------
for k = 1:nCases
    [ti, to] = ind2sub([n_tip_art, n_tip_art], total_ind_saved(k));
    t0 = tic;

    % Median-filter each frame (light smoothing)
    I = Image_all{k};
    I = arrayfun(@(z) medfilt2(I(:,:,z), [3,3]), 1:size(I,3), 'UniformOutput', false);
    I = cat(3, I{:});

    [H, W, T] = size(I);
    crop_sz   = 1;  % half-size around each centroid (results in up to 3x3 mean)

    % Clamp patch bounds
    r0 = max(round(centroid_shifted(:,2)) - crop_sz, 1);
    c0 = max(round(centroid_shifted(:,1)) - crop_sz, 1);
    r1 = min(round(centroid_shifted(:,2)) + crop_sz, H);
    c1 = min(round(centroid_shifted(:,1)) + crop_sz, W);

    % Mean over small patches per frame -> time-series per bead
    patch_means = cellfun(@(rs,re,cs,ce) ...
        mean(reshape(I(rs:re, cs:ce, :), [], T), 1), ...
        num2cell(r0), num2cell(r1), num2cell(c0), num2cell(c1), ...
        'UniformOutput', false);

    % [sub_num_beads x T]
    patch_means = cell2mat(patch_means);

    % Truncate/pad to N_REPEAT to be consistent
    if size(patch_means,2) >= N_REPEAT
        patch_means = patch_means(:, 1:N_REPEAT);
    else
        patch_means = [patch_means, zeros(num_beads, N_REPEAT - size(patch_means,2))];
    end

    % Save into the 4D cube
    Intensity_saved(:, ti, to, :) = patch_means;

    fprintf('Per-bead time-series: case %d/%d (ti=%d,to=%d) in %.2fs\n', ...
        k, nCases, ti, to, toc(t0));
end

% Save per-bead traces (optional diagnostic)
save(fullfile(save_folder, 'Intensity_saved_per_bead.mat'), 'Intensity_saved');

%% ============================ Patch Averaging =======================================
% Load patches saved previously (data_Folder/patches.mat)
loadpatchesPath = fullfile(data_Folder, 'patches.mat');
Intensity_saved = load(fullfile(save_folder, 'Intensity_saved_per_bead.mat')).Intensity_saved;

if ~exist(loadpatchesPath, 'file')
    error('Missing patches file: %s (run the first half to create it)', loadpatchesPath);
end
patches = load(loadpatchesPath).patches;

nPatches = numel(patches);
[~, nY, nX, ~] = size(Intensity_saved);
patchAvg = nan(nPatches, nY, nX, N_REPEAT);

for p = 1:nPatches
    idx = patches{p}.indices;     % bead indices in dim-1 of Intensity_saved
    if isempty(idx), continue; end
    % Mean across beads (dim 1) -> num_beads x nY x nX x N_REPEAT; squeeze -> nPatches x nY x nX x N_REPEAT
    patchAvg(p, :, :, :) = squeeze(mean(Intensity_saved(idx, :, :, :), 1));
end

% Drop the first few frames if needed 
Intensity_saved = patchAvg(:, :, :, 3:N_REPEAT);
clear patchAvg;  % free memory

% sub_num to "number of patches" from now on
sub_num = size(Intensity_saved, 1);

%% ============================ FFT & Response Maps ===================================
n_segment      = numel(target_segment);

% Frequency bins to sample
f = 20:1:100;
f = f(1:n_segment);

saved_fft_nonnormalized = zeros(sub_num, n_tip_art, n_tilt_art, n_segment);
saved_fft_normalized    = zeros(sub_num, n_tip_art, n_tilt_art, n_segment);
save_dc                 = zeros(n_tip_art, n_tilt_art);

for kk = 1:sub_num
    % nY x nX x T -> loop over cases
    Ipatch = squeeze(Intensity_saved(kk, :, :, :));  % [n_tip_art x n_tilt_art x T]
    for k = 1:nCases
        [ti, to] = ind2sub([n_tip_art, n_tip_art], total_ind_saved(k));

        % 1D FFT on the time trace at (ti,to)
        Y = abs(fft(squeeze(Ipatch(ti, to, :))));
        save_dc(ti, to) = Y(1);

        % Fill the maps at target frequencies (remap via tip/tilt index tables)
        for n_ind = 1:n_segment
            yi = f(n_ind) + 1;   % MATLAB 1-indexed FFT bin
            saved_fft_nonnormalized(kk, tip_ind_saved(n_ind, ti), tilt_ind_saved(n_ind, to), n_ind) = Y(yi);
            saved_fft_normalized(kk,    tip_ind_saved(n_ind, ti), tilt_ind_saved(n_ind, to), n_ind) = Y(yi) / max(Y(1), eps);
        end
    end
end


%% ============================ Peak Localization (Connected Components) ==============
% We detect bright blobs and use weighted centroids. If none found, fall back to
% normalized Gaussian fit.
% Locate peaks on each segment’s non-normalized FFT map and convert to angle.
pos_ref   = [ceil(n_tip_art/2), ceil(n_tilt_art/2)];

tip_tilt_saved = zeros(n_segment, 2);       % [dy, dx] in "pixels" then scaled by pix_angle
shift_saved    = zeros(n_segment, 2);

% Pixel->angle scale 
pix_angle = 0.15 / 0.6306;

% Precompute small grid for Gaussian fallback
[X, Y]  = meshgrid(-4:1:4, -4:1:4);
XY      = cat(3, X, Y);
gaussfn = @(B,XY) B(5) * exp(-((XY(:,:,1)-B(2))/B(3)).^2 - ((XY(:,:,2)-B(1))/B(4)).^2);

% Choose which "patch" to localize for (typical: best SNR). Here we localize for all patches,
for kk = 1:sub_num
    resnorm_save = [];

    for ii = 1:n_segment
        % Raw (non-normalized) map for stability
        M = squeeze(saved_fft_nonnormalized(kk, :, :, ii));
    
        % Normalize to [0, 1] for robust thresholding
        z = (M - min(M(:))) / max(max(M(:)) - min(M(:)), eps);

        % Threshold by mean + 1*std and keep top blobs by intensity
        thr  = mean(z(:)) + 1 * std(z(:));
        BW   = z > thr;
        BW   = bwareafilt(BW, [3, 40], 4);                           
        BW   = bwpropfilt(BW, z, 'MaxIntensity', 1, 'largest', 4);  
        props = regionprops(BW, z, 'Centroid', 'WeightedCentroid');

        % Fallback Gaussian fit if no blobs found
        if isempty(props)
            B0 = [0, 0, 0.5, 0.5, 1];   
            lb = [-4, -4, 0,   0,   0];
            ub = [ 4,  4, 4,   4,  10];
            try
                [B, resnorm] = lsqcurvefit(gaussfn, B0, XY, z, lb, ub);
            catch
                B       = [0,0,0,0,0];  
                resnorm = inf;
            end
            resnorm_save(end+1) = resnorm; 
            shift = B(1:2);  
        else
            % Use the brightest blob
            c = props(1).WeightedCentroid;     
            center = [ceil(n_tip_art/2), ceil(n_tilt_art/2)];
            shift  = c([2,1]) - center;        
        end

        % If fallback fit looks bad, discard
        if ~isempty(resnorm_save) && resnorm_save(end) > 4
            shift(:) = NaN;
        end

        shift_saved(ii, :)    = shift;
        tip_tilt_saved(ii, :) = shift * pix_angle;   

        if DEBUG_PLOT
            figure(10); subplot(4,5,ii); hold on;
            imagesc(z); axis image; colormap jet; colorbar; title(sprintf('Seg %d', ii)); 
            if all(isfinite(shift_saved(ii, :)))
                plot(shift_saved(ii,2)+pos_ref(2), shift_saved(ii,1)+pos_ref(1), ...
                     'kx', 'MarkerSize', 12, 'LineWidth', 2);
            end
            hold off;
        end
    end

    % Show progress
    fprintf('Progress: %d/%d (%.1f%%)\n', kk, sub_num, kk/sub_num*100);

    % Save per-patch result
    outName = sprintf('%s_centroid_patch_%03d.mat', ROUND_NAME, kk);
    % save(fullfile(save_folder, outName), 'tip_tilt_saved');
end
%%
function patches = divide_image_patches(I, centroids, cluster_dist, pad)
% DIVIDE_IMAGE_PATCHES  Split an image into patches based on clusters of bead centroids
%
%   patches = DIVIDE_IMAGE_PATCHES(I, centroids, cluster_dist, pad)
%
%   I             – input image (grayscale or RGB)
%   centroids     – M×2 array of [row, col] bead positions
%   cluster_dist  – maximum distance between beads to be considered same cluster (pixels)
%   pad           – number of pixels to pad around each cluster bounding box
%
%   patches       – cell array of structs with fields:
%                     .image – the image patch containing the cluster
%                     .bbox  – [r_min, r_max, c_min, c_max] in original image coords

    % Cluster centroids by proximity
    % Requires Statistics and Machine Learning Toolbox
    % Using single-linkage clustering with cutoff distance = cluster_dist
    labels = clusterdata(centroids, 'linkage', 'single', ...
                         'criterion', 'distance', 'cutoff', cluster_dist);
    nClusters = max(labels);

    % Initialize output
    patches = cell(nClusters,1);
    [H, W, ~] = size(I);

    for k = 1:nClusters
        % Find centroid indices in this cluster
        idx = find(labels == k);

        % Extract centroid positions
        rows = centroids(idx,1);
        cols = centroids(idx,2);

        % Compute bounding box with padding
        r_min = max(floor(min(rows)) - pad, 1);
        r_max = min(ceil(max(rows)) + pad, H);
        c_min = max(floor(min(cols)) - pad, 1);
        c_max = min(ceil(max(cols)) + pad, W);

        % Extract the patch
        patch_img = I(r_min:r_max, c_min:c_max, :);

        % Store results
        patches{k}.image   = patch_img;
        patches{k}.bbox    = [r_min, r_max, c_min, c_max];
        patches{k}.indices = idx;
    end
end
