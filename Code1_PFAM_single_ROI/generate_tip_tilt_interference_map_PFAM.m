%%
close all; clear; clc;
% Get current folder name
currentFolder = pwd;
[~, folderName] = fileparts(currentFolder);

% Check if folder name is "Code 1_FPAM"
if ~strcmp(folderName, 'Code1_PFAM_single_ROI')
    error('This script must be run from the "Code1_PFAM_single_ROI" folder. Current folder is: %s', folderName);
end

% Define paths
rootFolder = fileparts(currentFolder);
DM_data = fullfile(rootFolder, 'DM_data');
funcFolder = fullfile(rootFolder, 'function');
addpath(funcFolder);

% change folder to round 1, round 2 and round 3
save_folder = fullfile(currentFolder, 'data', 'round1');

% Load files
tip_ind_saved   = load(fullfile(save_folder, 'tip_ind_saved.mat')).tip_ind_saved;
tilt_ind_saved  = load(fullfile(save_folder, 'tilt_ind_saved.mat')).tilt_ind_saved;
total_ind_saved = load(fullfile(save_folder, 'total_ind_saved.mat')).total_ind_saved;

% Parameters
n_tip_art  = size(tip_ind_saved,2);
n_tilt_art = size(tilt_ind_saved,2);
nRepeat    = 200;


%% Show example image
Image_temp = double(loadtiff(fullfile(save_folder, 'Periodp16_ExposureTimep008_nFrames1_length_1000_7Segment_19_power_5_1\aberrated_Image.tif')));
Image_temp = Image_temp(2:size(Image_temp,1)-2,2:size(Image_temp,1)-2);
thresh = mean(Image_temp(:)) +5* std(Image_temp(:)); % use local threshold
Binary_Im = Image_temp > thresh;  %binary image
BW2 = bwareafilt(Binary_Im,5);

props = regionprops(BW2, Image_temp, 'Centroid', 'WeightedCentroid');
Binary_Im_new = zeros(size(Binary_Im,1),size(Binary_Im,2));
figure(10);
subplot(1,2,1); imagesc(Image_temp); axis image; colormap hot; colorbar; 
hold on
numObj = numel(props);
crop_temp =2;
for k = 1 : numObj
    plot(props(k).WeightedCentroid(1), props(k).WeightedCentroid(2), 'r*')
    plot(props(k).Centroid(1), props(k).Centroid(2), 'bo')
    Binary_Im_new(props(k).WeightedCentroid(2)-crop_temp:props(k).WeightedCentroid(2)+crop_temp,props(k).WeightedCentroid(1)-crop_temp:props(k).WeightedCentroid(1)+crop_temp)=1;
end

hold off;
Binary_Im_new = logical(Binary_Im_new);
subplot(1,2,2); imagesc(Binary_Im_new); colorbar; axis image; colormap hot;

Binary_Im_background = ~Binary_Im;

%% Setup
tip_tilt_fft = zeros(n_tip_art, n_tip_art);
save_value1  = zeros(n_tip_art, n_tip_art);
save_value2  = zeros(n_tip_art, n_tip_art);
save_value3  = zeros(n_tip_art, n_tip_art);
save_dc      = zeros(n_tip_art, n_tip_art);
Intensity_saved = zeros(n_tip_art, n_tilt_art, nRepeat);   % shape matches your earlier code

% Common path prefix for each (tip,tilt) case
save_name_full = strcat(save_folder, ...
    '\Periodp16_ExposureTimep008_nFrames1_length_1000_7Segment_19_power_5_');

%% --------- Preload all image stacks and references to avoid repeated I/O ---------
nCases = length(total_ind_saved);
Image_ref_all = cell(nCases, 1);
Image_all     = cell(nCases, 1);

for ii = 1:nCases
    [tip_idx, tilt_idx] = ind2sub([n_tip_art, n_tip_art], total_ind_saved(ii));
    save_name  = strcat(save_name_full, num2str((tip_idx-1)*n_tip_art + tilt_idx));
    final_path = fullfile(save_name);

    % Load once per case
    Image_ref_all{ii} = double(loadtiff(fullfile(final_path, 'aberrated_Image.tif')));
    Image_all{ii}     = double(loadtiff(fullfile(final_path, 'test.tiff')));
end

%% --------- Process each case (vectorized over frames) ---------
for ii = 1:nCases
    [n_tip_art_ind, n_tilt_art_ind] = ind2sub([n_tip_art, n_tip_art], total_ind_saved(ii));

    t0 = tic;

    % Pull preloaded stacks
    Image_ref = Image_ref_all{ii};
    Image     = Image_all{ii};

    % (Optional) crop borders as in your original code
    % Keep single z for ref; keep all frames for Image
    Image_ref = Image_ref(2:end-2, 2:end-2, 1);
    Image     = Image(    2:end-2, 2:end-2, :);

    % ----- Build ROI mask (Binary_Im_new) from Image_ref once -----
    % Threshold + keep largest N blobs
    thresh   = mean(Image_ref(:)) + 5*std(Image_ref(:));
    BW       = Image_ref > thresh;
    BW2      = bwareafilt(BW, 5);  % keep 5 objects

    % Weighted centroids for 5 objects
    props    = regionprops(BW2, Image_ref, 'WeightedCentroid');

    % Draw 5×5 (crop_sz=2) squares centered at the (floored) weighted centroid positions
    [H, W] = size(Image_ref);
    crop_sz = 2;  % 5x5 region
    ROI = false(H, W);
    for k = 1:numel(props)
        r0 = max(1, floor(props(k).WeightedCentroid(2)) - crop_sz);
        c0 = max(1, floor(props(k).WeightedCentroid(1)) - crop_sz);
        r1 = min(H, floor(props(k).WeightedCentroid(2)) + crop_sz);
        c1 = min(W, floor(props(k).WeightedCentroid(1)) + crop_sz);
        ROI(r0:r1, c0:c1) = true;
    end
    ROI_bg = ~ROI;

    % ----- Vectorized per-frame metric over all frames -----
    T = size(Image, 3);
    roi_count = nnz(ROI);
    bg_count  = nnz(ROI_bg);

    % Sum over x,y for every frame in one pass
    roi_sum = squeeze(sum(sum(Image .* ROI,   1), 2));   % [T x 1]
    bg_sum  = squeeze(sum(sum(Image .* ROI_bg,1), 2));   % [T x 1]

    Metric     = roi_sum ./ max(roi_count, 1);   % avoid divide by zero
    Background = bg_sum  ./ max(bg_count,  1);

    % Choose which metric to save (your original used ROI mean only)
    % Intensity = Metric - Background;  % contrast metric (optional)
    Intensity = Metric;

    % Ensure length compatibility with nRepeat (truncate/pad if needed)
    if numel(Intensity) < nRepeat
        Intensity = [Intensity(:); zeros(nRepeat - numel(Intensity), 1)];
    elseif numel(Intensity) > nRepeat
        Intensity = Intensity(1:nRepeat);
    end

    % Save
    Intensity_saved(n_tip_art_ind, n_tilt_art_ind, :) = Intensity;

    % Diagnostics
    dt = toc(t0);
    fprintf('(%d/%d) tip=%d, tilt=%d, dt=%.3fs\n', ii, nCases, n_tip_art_ind, n_tilt_art_ind, dt);
end

% Save once at the end
save(fullfile(save_folder, "Intensity_saved.mat"), "Intensity_saved");

%%
Intensity_saved = load(fullfile(save_folder, '\Intensity_saved.mat')).Intensity_saved;
% Intensity_saved = load(fullfile(save_name_temp, '\Intensity_saved_new_aberrated.mat')).Intensity_saved;
Intensity_saved = Intensity_saved(:,:,3:nRepeat);

target_segment =[1 20 23 26 29 32 35 92 95 98 101 104 107 110 113 116 119 122 125]; 
% target_segment =[9 13 17 44 52 60 63 66 73 76 83 86 137 140 151 154 165 168];
% target_segment = [11 15 19 40 48 56 68 71 78 81 88 91 130 133 144 147 158 161];
f = [20:1:100]'
f = f(1:length(target_segment))

n_segment = length(target_segment);
start_fs = 100; step = 4;
saved_fft_nonnormalized = zeros(n_tip_art, n_tilt_art, n_segment);
saved_fft_normalized = zeros(n_tip_art, n_tilt_art, n_segment);

% for n_tip_art_ind = 1:n_tip_art
% for n_tilt_art_ind = 1:n_tilt_art
for n_total_art_ind = 1:length(total_ind_saved)
% for n_total_art_ind = 40

[n_tip_art_ind, n_tilt_art_ind] = ind2sub([n_tip_art, n_tilt_art], total_ind_saved(n_total_art_ind));

% for n_tip_art_ind = 3
% for n_tilt_art_ind = 9
disp([tip_ind_saved(1, n_tip_art_ind), tilt_ind_saved(1, n_tilt_art_ind)]);
% fft
Y = abs(fft(squeeze(Intensity_saved(n_tip_art_ind, n_tilt_art_ind, :))));
Y_show = Y; 
Y_show(1:5)=0;
% Y(1:10)=0;
figure(1122);
plot(squeeze(Intensity_saved(n_tip_art_ind, n_tilt_art_ind, :)));
figure(1123);
plot(Y_show);
for n_ind= 1:n_segment
    saved_fft_normalized(tip_ind_saved(n_ind, n_tip_art_ind), tilt_ind_saved(n_ind, n_tilt_art_ind), n_ind) = Y(f(n_ind)+1)./Y(1);
    saved_fft_nonnormalized(tip_ind_saved(n_ind, n_tip_art_ind), tilt_ind_saved(n_ind, n_tilt_art_ind), n_ind) = Y(f(n_ind)+1);


end

save_dc(n_tip_art_ind,n_tilt_art_ind) = Y(1);
end
%%
figure(1124);
plot(Y_show, '-o'); xlim([21,32]); xticks([21:2:32])

%%
figure(12311);
for ii=1:length(target_segment)
subplot(4,5,ii); imagesc(medfilt2(saved_fft_normalized(:, :, ii))); axis image;  colormap gray; title([ii]); colorbar;
end
figure(12312);
for ii=1:length(target_segment)
subplot(4,5,ii); imagesc(medfilt2(saved_fft_nonnormalized(:, :, ii))); axis image;  colormap gray; title([ii]); colorbar;
end
figure(123133);
for ii=1:length(target_segment)
subplot(4,5,ii); imagesc((saved_fft_nonnormalized(:, :, ii))); axis image;  colormap jet; title([ii]); colorbar;
end

figure(12313);
imagesc(save_dc); axis image;  colormap gray; title([ii]); colorbar;
save(strcat(save_folder, '\tip_tilt_map.mat'), 'saved_fft_nonnormalized');
%% find peak
test_centroid = medfilt2(saved_fft_normalized(:, :, 1), [3,3]);
test_centroid = test_centroid - min(test_centroid(:));
test_centroid = test_centroid ./ max(test_centroid(:)) *255;

test_centroid_label = bwlabel(test_centroid);
figure(1231);
subplot(1,2,1); imagesc(test_centroid); axis image; colormap gray;
subplot(1,2,2); imagesc(test_centroid_label); axis image; colormap gray;


% props = regionprops(test_centroid_label, test_centroid, 'Centroid', 'WeightedCentroid')
% [6 6]- props.WeightedCentroid
pos_ref=[ceil(n_tip_art/2),ceil(n_tilt_art/2)];
pos0=[ceil(n_tip_art/2),ceil(n_tilt_art/2)];
Period = n_tip_art;
isDisplay = true;
sig = 1;
ItN =100;
tip_tilt_saved = zeros(length(target_segment),2);
shift_saved = zeros(length(target_segment),2);

[X, Y] = meshgrid(-floor(Period/2):1:floor(Period/2), -floor(Period/2):1:floor(Period/2));
XY(:,:,1) = X;
XY(:,:,2) = Y;

gaussfit = @(B,XY)  B(5) * exp(-((XY(:,:,1) - B(2))/B(3)).^2 - ((XY(:,:,2)-B(1))/B(4)).^2);
Xdata = [X(:), Y(:)];

pix_angle = 0.15 ./ 0.6306;

for ii=1:length(target_segment)
% 
    test_data_med = saved_fft_nonnormalized(:, :, ii);
    z = (test_data_med - min(test_data_med(:))) / (max(test_data_med(:)) - min(test_data_med(:)));
    thresh = mean(z(:))+1*std(z(:))
    test_data_label =  z > thresh;
    BW2 = bwareafilt(test_data_label,[1, 80],4);

    BW2 = bwpropfilt(BW2,z,'MaxIntensity',1, "largest",4);

    props = regionprops(BW2, z, 'Centroid', 'WeightedCentroid');
    % tip/tilt interference map quality check using gaussian fit
    if length(props)==0
        B = lsqcurvefit(gaussfit, [0, 0, 0.5, 0.5, 0], XY, z, [-floor(Period/2), -floor(Period/2), 0, 0, 0], [floor(Period/2), floor(Period/2), floor(Period/2), floor(Period/2), 10]);
        % shift_saved(ii,:) = B(1:2);
        shift_saved(ii,:) = [NaN, NaN];

    else
        shift_temp = props(1).WeightedCentroid;
        shift_saved(ii,:) = shift_temp([2,1])-[ceil(n_tip_art/2),ceil(n_tilt_art/2)];

    end
        
    figure(12314);
    subplot(4,5,ii); imagesc(z); axis image;  colormap jet; title([ii]); colorbar; hold on; %caxis([0,0.025]);
    scatter(shift_saved(ii,2)+pos_ref(1),shift_saved(ii,1)+pos_ref(2), 500, 'x', 'Color', 'g');
    hold off
    tip_tilt_saved(ii,:) = shift_saved(ii,:)*pix_angle;
    figure(12315);
    subplot(4,5,ii); imagesc(BW2); axis image;  colormap gray; title([num2str(ii), ' ',num2str(thresh)]); colorbar; hold on;
    hold off
    figure(12316);
    scatter(shift_saved(ii,2)+pos_ref(1),shift_saved(ii,1)+pos_ref(2), 100, 'x', 'Color', 'b'); axis([0 n_tip_art 0 n_tilt_art]); hold on; 
end
%%
% save(strcat(save_name_temp, '\round1_1_downsampled8.mat'), 'tip_tilt_saved')
% save(strcat(save_name_temp, '\round1_2_quantile_p5.mat'), 'tip_tilt_saved')
save(strcat(save_folder, '\round3_2_1p0_noft_crop2_single.mat'), 'tip_tilt_saved')
% save(strcat(save_name_temp, '\round2_2_p0.mat'), 'tip_tilt_saved')


% save(strcat(save_name_temp, '\round1_2_square_centroid.mat'), 'tip_tilt_saved')
% save(strcat(save_name_temp, '\round2_2_gaussian.mat'), 'tip_tilt_saved')
% % 
%%
% pos_ref=[ceil(n_tip_art/2) ceil(n_tilt_art/2)];
% pos0=[ceil(n_tip_art/2) ceil(n_tilt_art/2)];
% Period = n_tip_art;
% isDisplay = true;
% sig = 1;
% ItN =100;
% tip_tilt_saved = zeros(length(target_segment),2);
% shift_saved = zeros(length(target_segment),2);
% 
% pix_angle = 0.12 ./ 0.6306;
% % pix_angle = 0.24 ./ 0.6306;
% 
% for ii=1:length(target_segment)
% % for ii=1:20
% % for ii=21:length(target_segment)
%     test_data_med = padarray(saved_fft_normalized(:, :, ii), [3,3],'both');
% 
%     % test_data_med = medfilt2(saved_fft_normalized(:, :, ii), [3,3]);
%     test_data_med = imgaussfilt(test_data_med, 0.5);
% 
%     % test_data_med = saved_fft_normalized(:, :, ii);
%     pos_ref=[ceil(size(test_data_med,1)/2) ceil(size(test_data_med,2)/2)];
%     pos0=[ceil(size(test_data_med,1)/2) ceil(size(test_data_med,2)/2)];
%     Period = size(test_data_med,1);
% 
%     % if ii>20
%     %     test_data_med = medfilt2(saved_fft_normalized(:, :, ii), [3,3]);
%     % end
%     thresh = mean(test_data_med(:))+ 1*std(test_data_med(:));
%     % test_data_label =  test_data_med > thresh;
%     % test_data_med(~test_data_label)=0;
%     % test_data_med = test_data_med - min(test_data_med(:));
%     % test_data_med = test_data_med ./ max(test_data_med(:)) *255;
%     shift = single_shift_finder(test_data_med, pos0, pos_ref, Period, isDisplay, sig, ItN);
% 
%     shift_saved(ii,:) = shift([2,1]);
%     figure(12314);
%     subplot(6,6,ii); imagesc(test_data_med); axis image;  colormap gray; title([ii]); colorbar; hold on;
%     scatter(shift_saved(ii,2)+pos_ref(1),shift_saved(ii,1)+pos_ref(2), 500, 'x', 'Color', 'g');axis([0 n_tip_art 0 n_tilt_art]);
%     hold off
%     tip_tilt_saved(ii,:) = shift([2,1])*pix_angle;
% 
% end
% %%
% tip_tilt_saved_f_20 = load('D:\WidefieldAO\240515\tip_tilt_forloop26\tip_tilt_saved_f_20.mat').tip_tilt_saved;
% tip_tilt_saved_f_40 = load('D:\WidefieldAO\240515\tip_tilt_forloop26\tip_tilt_saved_f_40.mat').tip_tilt_saved;
% 
% tip_tilt_saved = tip_tilt_saved_f_20;
% tip_tilt_saved([1:7, 8 11 14 17],:) = tip_tilt_saved_f_40([1:7, 8 11 14 17],:);
% %%
% make_zeros = [8 9];
% tip_tilt_saved(make_zeros, :)=repmat([NaN,NaN], length(make_zeros),1)
% % tip_tilt_saved(make_zeros, :)=repmat([0,0], length(make_zeros),1)
% 
% %%
% tip_tilt_saved_save_center =  mean(tip_tilt_saved, 1, 'omitnan');
% temp = tip_tilt_saved - tip_tilt_saved_save_center;
% temp = temp*0.6306;
% error = temp + zeros(size(temp,1),2);
% dist_mean = mean(sqrt(sum(error.^2,2)))
% 
% %%
% close all; clear; clc;
% Param.level = 7; % 7
% Param.SegN = Param.level*(Param.level+1)*3+1;
% Param.UnreachableSeg = [128, 135, 142, 149, 156, 163]; 
% Param.LockedSeg = [5, 70,  80, 100, 155]; % 20190924: added #80 
% 
% Param.Outlier_hex = [169, 128, 129,  134, 135, 136,  141, 142, 143, ...
%                   148, 149, 150,  155, 156, 157, 162, 163, 164, ...
%                  5, 70, 80, 100]; 
% 
% focal_len = 46;        % focal-length of lenslet-array, milimeter
% SHcam_pitch = 6.5;     % pitch size of SH-camera, micron
% 
% DM_pitch = 606;    % pitch size of SM, micron
% Mag_DM2LA = 85/300*175/60;      % ~0.83
% SH_ref_name = 'SHref-10162023'
% load(['data\' SH_ref_name '.mat']);                            % load the SH shift reference: PosRef
% load('data\IndexTransfer.mat', 'IndexSeq_hex2jk', 'IndexCell_jk2hex');       % index convertion
% Convert_hex2jk = IndexSeq_hex2jk;
% Convert_jk2hex = IndexCell_jk2hex;
% SH_name = 'SH1';
% SH_im = double(imread(['data\' SH_name '.tif']));              % load SH image
% 
% % load(['data\PR3_240308.mat']);
% % IrisAO_PTT0_1 = IrisAO_PTT;
% % IrisAO_PTT0_2 = IrisAO_PTT;
% 
% load(['data\IrisAO_PTT_2_DM_wrap_0625_agarose1_2hexagonal_iter1.mat']);
% IrisAO_PTT0_1 = IrisAO_PTT_2_DM_wrap;
% IrisAO_PTT0_2 = IrisAO_PTT_2_DM_wrap;
% % % 
% IrisAO_PTT0_1 = zeros(169,3);
% IrisAO_PTT0_2 = zeros(169,3);
% 
% % IrisAO_PTT_0 = zeros(169,3);
% % load(['data\PR3_0_768_240312.mat']);
% % load(['data\pr3_0_1024_240618_2.mat']);
% % load(['data\IrisAO_PTT_2_DM_wrap_0523_artificial_rand_p02.mat']);
% % load(['data\IrisAO_PTT_2_DM_wrap_0523_artificial_asti_p1.mat']);
% % load(['data\IrisAO_PTT_2_DM_wrap_0604_syscorr_3hexagonal_iter2.mat']);
% % load(['data\pr1_1900_1024.mat']);
% % load(['data\PR3_240425.mat']);
% load(['data\pr3_agarose1_240625.mat']);
% 
% % pr1_0_1024_240611
% % IrisAO_PTT_GT =IrisAO_PTT_2_DM_wrap;
% IrisAO_PTT_GT =IrisAO_PTT;
% % IrisAO_PTT_GT =IrisAO_PTT_2_DM_wrap;
% 
% load(['data\ZernikeCoeff2PTT2.mat']);
% IrisAO_PTT_GT =  -0.05*squeeze(PTT_mode(4,:,:));
% 
% 
% 
% %%
% tip_tilt_total1 = NaN(169,2);
% tip_tilt_total2 = NaN(169,2);
% tip_tilt_total3 = NaN(169,2);
% 
% round1_segment = [];
% round2_segment = [];
% round3_segment = [];
% 
% round1_tip_tilt = load("D:\WidefieldAO\240701\tip_tilt_forloop5\tip_tilt_saved_round1.mat").tip_tilt_saved;
% round2_tip_tilt = load("D:\WidefieldAO\240701\tip_tilt_forloop8\tip_tilt_saved_round2.mat").tip_tilt_saved;
% round3_tip_tilt = load("D:\WidefieldAO\240701\tip_tilt_forloop8\tip_tilt_saved_round2.mat").tip_tilt_saved;
% 
% 
% % round12_tip_tilt = [round1_tip_tilt; round2_tip_tilt];
% % % round12_tip_tilt =  load("D:\WidefieldAO\240416\tip_tilt_forloop1\tip_tilt_saved_round1.mat").tip_tilt_saved;
% % % 
% % round12_centered = mean(round12_tip_tilt, 1, 'omitnan');
% % round1_tip_tilt = round1_tip_tilt - round12_centered;
% % round2_tip_tilt = round2_tip_tilt - round12_centered;
% % round1_tip_tilt = round1_tip_tilt- round12_centered;
% 
% % round123_tip_tilt = [round1_tip_tilt;];
% round123_tip_tilt = [round1_tip_tilt; round2_tip_tilt;];
% % round123_tip_tilt = [round1_tip_tilt; round2_tip_tilt; round3_tip_tilt];
% round123_centered = mean(round123_tip_tilt, 1, 'omitnan');
% % round123_centered=0;
% round1_tip_tilt = round1_tip_tilt - round123_centered;
% round2_tip_tilt = round2_tip_tilt - round123_centered;
% round3_tip_tilt = round3_tip_tilt - round123_centered;
% % round1_tip_tilt = [-0.0032 0.0239]./0.6303;
% % round123_tip_tilt = load("D:\WidefieldAO\240416\tip_tilt_forloop1\tip_tilt_saved_round1.mat").tip_tilt_saved;
% % round123_centered = mean(round123_tip_tilt, 1, 'omitnan');
% 
% % round1_tip_tilt = round1_tip_tilt ;
% % round2_tip_tilt = round2_tip_tilt ;
% % round3_tip_tilt = round3_tip_tilt ;
% 
% % target_segment1 =[1 20 23 26 29 32 35 92 95 98 101 104 107 110 113 116 119 122 125]; 
% % target_segment2 =[9 13 17 44 52 60 63 66 73 76 83 86 137 140 151 154 165 168];
% % target_segment3 = [11 15 19 40 48 56 68 71 78 81 88 91 130 133 144 147 158 161];
% % target_segment1 =[2 14 38 43 46 54 57 77 96 124 145 153];
% % target_segment2 =[11 17 40 48 52 60 98 101 104 116 119 122 130 147 151 168];
% 
% % target_segment2 =[11 17 41 47 53 59 94 98 101 104 108 112 116 119 122 126];
% 
% target_segment1 = [1 20 23 26 29 32 35];
% % target_segment2 = [62 95 67 101 72 107 77 113 82 119 87 125];
% % target_segment3 = [62 95 67 101 72 107 77 113 82 119 87 125];
% target_segment2 = [92 95 98 101 104 107 110 113 116 119 122 125];
% target_segment3 = [92 95 98 101 104 107 110 113 116 119 122 125];
% 
% % target_segment1 =[9 13 17 44 52 60];
% % target_segment2 =[63 66 73 76 83 86 137 140 151 154 165 168];
% % target_segment3 =[63 66 73 76 83 86 137 140 151 154 165 168];
% 
% % target_segment1 =[1 20 23 26 29 32 35 92 95 98 101 104 107 110 113 116 119 122 125]; 
% % target_segment2 =[1 20 23 26 29 32 35 92 95 98 101 104 107 110 113 116 119 122 125]; 
% % target_segment3 =[1 20 23 26 29 32 35 92 95 98 101 104 107 110 113 116 119 122 125]; 
% 
% % target_segment1 =[1 20 23 26 29 32 35 92 95 98 101 104 107 110 113 116 119 122 125]; 
% % target_segment2 =[9 13 17 44 52 60 63 66 73 76 83 86 137 140 151 154 165 168]; 
% % target_segment3 =[11 15 19 40 48 56 68 71 78 81 88 91 130 133 144 147 158 161]; 
% 
% % target_segment1 = [92 95 98 101 104 107 110 113 116 119 122 125];
% % target_segment2 = [44 52 60 63 66 73 76 83 86];
% % target_segment3 = [11 15 19 40 48 56];
% for ii=1:length(target_segment1)
%     % disp(ii)
%     [neighbours_hex, ~] = Find6NeighbrIdx(target_segment1(ii), Param, Convert_hex2jk, Convert_jk2hex);
%     if target_segment1(ii) ==1 || target_segment1(ii) == 2
%         round1_segment = [round1_segment; [target_segment1(ii); neighbours_hex]'];
%         % [ii; neighbours_hex]'
%         tip_tilt_total1([target_segment1(ii); neighbours_hex],:) = repmat(round1_tip_tilt(ii,:),7,1);
%     else
%         neighbours_hex(neighbours_hex==1) = 163;
%         round1_segment = [round1_segment; [target_segment1(ii); neighbours_hex]'];
%         % [ii; neighbours_hex]'
%         tip_tilt_total1([target_segment1(ii); neighbours_hex],:) = repmat(round1_tip_tilt(ii,:),7,1);
%     end
% end
% for ii=1:length(target_segment2)
%     % disp(ii)
%     [neighbours_hex, ~] = Find6NeighbrIdx(target_segment2(ii), Param, Convert_hex2jk, Convert_jk2hex);
%     if target_segment2(ii) ==1 || target_segment2(ii) == 2
%         round2_segment = [round2_segment; [target_segment2(ii); neighbours_hex]'];
%         % [ii; neighbours_hex]'
%         tip_tilt_total2([target_segment2(ii); neighbours_hex],:) = repmat(round2_tip_tilt(ii,:),7,1);
%     else
%         neighbours_hex(neighbours_hex==1) = 163;
%         round2_segment = [round2_segment; [target_segment2(ii); neighbours_hex]'];
%         % [ii; neighbours_hex]'
%         tip_tilt_total2([target_segment2(ii); neighbours_hex],:) = repmat(round2_tip_tilt(ii,:),7,1);
%     end
% end
% 
% for ii=1:length(target_segment3)
%     % disp(ii)
%     [neighbours_hex, ~] = Find6NeighbrIdx(target_segment3(ii), Param, Convert_hex2jk, Convert_jk2hex);
%     if target_segment3(ii) ==1 || target_segment3(ii) == 2
%         round3_segment = [round3_segment; [target_segment3(ii); neighbours_hex]'];
%         % [ii; neighbours_hex]'
%         tip_tilt_total3([target_segment3(ii); neighbours_hex],:) = repmat(round3_tip_tilt(ii,:),7,1);
%     else
%         neighbours_hex(neighbours_hex==1) = 163;
%         round3_segment = [round3_segment; [target_segment3(ii); neighbours_hex]'];
%         % [ii; neighbours_hex]'
%         tip_tilt_total3([target_segment3(ii); neighbours_hex],:) = repmat(round3_tip_tilt(ii,:),7,1);
%     end
% end
% % for ii=1:length(target_segment3)
% %     % disp(ii)
% %     [neighbours_hex, ~] = Find6NeighbrIdx(target_segment3(ii), Param, Convert_hex2jk, Convert_jk2hex);
% %     neighbours_hex(neighbours_hex==1) = 163;
% %     used_segment_total = [used_segment_total; [target_segment3(ii); neighbours_hex]'];
% %     % [ii; neighbours_hex]'
% %     tip_tilt_total3([target_segment3(ii); neighbours_hex],:) = repmat(round3_tip_tilt(ii,:),7,1);
% % end
% 
% 
% % used_segment_total = [round1_segment;];
% used_segment_total = [round1_segment; round2_segment];
% % used_segment_total = [round1_segment; round2_segment; round3_segment];
% 
% used_segment_total = unique(used_segment_total(:));
% %%
% filepath = 'D:\WidefieldAO\data\';
% load([filepath 'ZernikeCoeff2PTT.mat'], 'PTT_mode');
% PTT_mode_round1 = PTT_mode;
% PTT_mode_round2 = PTT_mode;
% PTT_mode_round3 = PTT_mode;
% 
% for temp=1:length(target_segment1)
% 
%     PTT_mode_round1(2,round1_segment(temp,:),:) = [0.00204724857001837,-1.20486197704087e-18,0.630552559565678;...
%     0.00204724857001847,-4.53257600886802e-17,0.630552559565683;...
%     0.276378556952489,9.79897033613380e-16,0.630552559565678;...
%     0.276378556952489,-9.47767380892290e-16,0.630552559565678;...
%     0.00204724857001847,0,0.630552559565678;...
%     -0.272284059812451,9.99232199625893e-16,0.630552559565673;...
%     -0.272284059812451,-1.00336315497575e-15,0.630552559565672];
% % 
%     PTT_mode_round1(3,round1_segment(temp,:),:) = [0.00204724857001844,0.630552559565684,2.26820881735948e-17;...
%     0.317323528352860,0.630552559565684,9.63736666622372e-16;...
%     0.157638139891421,0.630552559565682,9.88959878367924e-16;...
%     -0.157638139891421,0.630552559565682,-6.92432182390090e-16;...
%     -0.317323528352860,0.630552559565682,0;...
%     -0.157638139891421,0.630552559565682,-9.03454401488104e-17;...
%     0.157638139891421,0.630552559565682,4.12648576278340e-16];
% 
% end
% 
% for temp=1:length(target_segment2)
% 
%     PTT_mode_round2(2,round2_segment(temp,:),:) = [0.00204724857001837,-1.20486197704087e-18,0.630552559565678;...
%     0.00204724857001847,-4.53257600886802e-17,0.630552559565683;...
%     0.276378556952489,9.79897033613380e-16,0.630552559565678;...
%     0.276378556952489,-9.47767380892290e-16,0.630552559565678;...
%     0.00204724857001847,0,0.630552559565678;...
%     -0.272284059812451,9.99232199625893e-16,0.630552559565673;...
%     -0.272284059812451,-1.00336315497575e-15,0.630552559565672];
% % 
%     PTT_mode_round2(3,round2_segment(temp,:),:) = [0.00204724857001844,0.630552559565684,2.26820881735948e-17;...
%     0.317323528352860,0.630552559565684,9.63736666622372e-16;...
%     0.157638139891421,0.630552559565682,9.88959878367924e-16;...
%     -0.157638139891421,0.630552559565682,-6.92432182390090e-16;...
%     -0.317323528352860,0.630552559565682,0;...
%     -0.157638139891421,0.630552559565682,-9.03454401488104e-17;...
%     0.157638139891421,0.630552559565682,4.12648576278340e-16];
% 
% end
% 
% for temp=1:length(target_segment3)
% 
%     PTT_mode_round2(2,round3_segment(temp,:),:) = [0.00204724857001837,-1.20486197704087e-18,0.630552559565678;...
%     0.00204724857001847,-4.53257600886802e-17,0.630552559565683;...
%     0.276378556952489,9.79897033613380e-16,0.630552559565678;...
%     0.276378556952489,-9.47767380892290e-16,0.630552559565678;...
%     0.00204724857001847,0,0.630552559565678;...
%     -0.272284059812451,9.99232199625893e-16,0.630552559565673;...
%     -0.272284059812451,-1.00336315497575e-15,0.630552559565672];
% % 
%     PTT_mode_round2(3,round3_segment(temp,:),:) = [0.00204724857001844,0.630552559565684,2.26820881735948e-17;...
%     0.317323528352860,0.630552559565684,9.63736666622372e-16;...
%     0.157638139891421,0.630552559565682,9.88959878367924e-16;...
%     -0.157638139891421,0.630552559565682,-6.92432182390090e-16;...
%     -0.317323528352860,0.630552559565682,0;...
%     -0.157638139891421,0.630552559565682,-9.03454401488104e-17;...
%     0.157638139891421,0.630552559565682,4.12648576278340e-16];
% 
% end
% 
% round1_PTT = tip_tilt_total1(:,1).*squeeze(PTT_mode_round1(2,:,:))+tip_tilt_total1(:,2).*squeeze(PTT_mode_round1(3,:,:));
% round2_PTT = tip_tilt_total2(:,1).*squeeze(PTT_mode_round2(2,:,:))+tip_tilt_total2(:,2).*squeeze(PTT_mode_round2(3,:,:));
% round3_PTT = tip_tilt_total3(:,1).*squeeze(PTT_mode_round3(2,:,:))+tip_tilt_total3(:,2).*squeeze(PTT_mode_round3(3,:,:));
% 
% %%
% % PTT_total = zeros(1,169,3);
% PTT_total = zeros(2,169,3);
% % PTT_total = zeros(3,169,3);
% 
% PTT_total(1,:,:)= round1_PTT;
% PTT_total(2,:,:)= round2_PTT;
% % PTT_total(3,:,:)= round3_PTT;
% 
% % tip_tilt_total(3,:,:)= tip_tilt_total3;
% 
% PTT_total_mean = squeeze(mean(PTT_total, 1, "omitnan"));
% % PTT_total_mean = IrisAO_PTT;
% 
% %
% % tip_tilt_total = zeros(2,169,2);
% % tip_tilt_total(1,:,:)= tip_tilt_total1;
% % tip_tilt_total(2,:,:)= tip_tilt_total2;
% % % tip_tilt_total(3,:,:)= tip_tilt_total3;
% % 
% % tip_tilt_total_mean = squeeze(mean(tip_tilt_total, 1, "omitnan"));
% % %%
% % tip_tilt_total_mean = 0.05*squeeze(PTT_mode(8,:,3:-1:2));
% % % tip_tilt_total_mean = -0.4*ones(169,2);
% % % tip_tilt_total
% % %%
% % slopes = zeros(169,2);
% % slopes(:,1)= -2/1000 *tip_tilt_total_mean(:, 1);
% % slopes(:,2)= -2/1000 *tip_tilt_total_mean(:, 2);
% % slopes_new = slopes;
% % unused_segment_total = setdiff([1:169], used_segment_total(:));
% % slopes_new(unused_segment_total,:) = 0;
% % flag = find(isnan(slopes_new(:,1)));
% % 
% % % disp(Param.Outlier_hex)
% % flag = find(isnan(slopes_new(:,1)));
% % slopes_new(flag,:,:)=0;
% % unused_segment_total = [unused_segment_total, flag'];
% % 
% % slopes_new = CalcLockedSegShift(slopes_new, Param.Outlier_hex', Param, Convert_hex2jk, Convert_jk2hex);
% 
% 
% 
% %
% flag_nan = find(isnan(PTT_total_mean(:,1)));
% 
% Param.Outlier_hex = union(Param.Outlier_hex,flag_nan);
% disp(Param.Outlier_hex)
% %
% PTT_total_mean(flag_nan,:,:)=0;
% unused_segment_total = setdiff([1:169], used_segment_total(:));
% 
% unused_segment_total = [unused_segment_total, flag_nan'];
% 
% PTT_total_mean = CalcLockedSegShift(PTT_total_mean, unused_segment_total, Param, Convert_hex2jk, Convert_jk2hex);
% Seg_active = 1:Param.SegN;
% PTT_total_mean(Seg_active, 1) = PTT_total_mean(Seg_active, 1) - PTT_total_mean(1, 1);
% PTT_total_mean(Seg_active, 2) = PTT_total_mean(Seg_active, 2) - PTT_total_mean(1, 2);
% PTT_total_mean(Seg_active, 3) = PTT_total_mean(Seg_active, 3) - PTT_total_mean(1, 3);
% %
% slopes(:,1)= -2/1000 *PTT_total_mean(:, 3);
% slopes(:,2)= -2/1000 *PTT_total_mean(:, 2);
% slopes_new = slopes;
% 
% 
% %
% slopes_new = slopes;
% flag_nan = find(isnan(slopes_new(:,1)));
% 
% Param.Outlier_hex = union(Param.Outlier_hex,flag_nan);
% % disp(Param.Outlier_hex)
% slopes_new(flag_nan,:)=0;
% unused_segment_total = [unused_segment_total, flag_nan'];
% 
% % slopes_new = CalcLockedSegShift(slopes_new, Param.Outlier_hex', Param, Convert_hex2jk, Convert_jk2hex);
% % Seg_active = 1:Param.SegN;
% % slopes_new(Seg_active, 1) = slopes_new(Seg_active, 1) - slopes_new(1, 1);
% % slopes_new(Seg_active, 2) = slopes_new(Seg_active, 2) - slopes_new(1, 2);
% %
% % slopes_new([135, 136, 137, 138],:) = repmat(slopes_new(99,:),4,1);
% % slopes_new([138,101, 70],:) = repmat(slopes_new(45,:),3,1);
% % slopes_new([137,100, 69],:) = repmat(slopes_new(44,:),3,1);
% % slopes_new([136],:) = repmat(slopes_new(99,:),1,1);
% % slopes_new([135 98],:) = repmat(slopes_new(67,:),2,1);
% 
% % slopes_new([100, 101],:) = repmat(slopes_new(69,:),2,1);
% % slopes_new(98,:) = slopes_new(97,:);
% % slopes_new = CalcLockedSegShift(slopes_new, flag', Param, Convert_hex2jk, Convert_jk2hex);
% 
% %%
% slopes = slopes_new;
% W = zeros(Param.SegN, 1);   %  current wavefront, initialized to zero
% W_update = zeros(Param.SegN, 1);    %  updated wavefront
% 
% % W = PTT_total_mean(:, 1);   %  current wavefront, initialized to zero
% % W_update = PTT_total_mean(:, 1);    %  updated wavefront
% % W_update = zeros(Param.SegN, 1);    %  updated wavefront
% 
% tol = 1e-3;     % tolerance for convergency % default: 1e-3
% iterN = 0;
% delta_W = [];
% Seg_active = 1:Param.SegN;
% % 20210722, for letting the reconstruction converge for fewer segments
% % we also exclude outliers for phase reconstruction
% Param.Outlier_hex = [169, 128, 129,  134, 135, 136,  141, 142, 143, ...
%                   148, 149, 150,  155, 156, 157, 162, 163, 164, ...
%                  5, 70, 80, 100]; 
% 
% Seg_active(flag_nan) = [];   % Should we remove them?? --> No
% % disp(Seg_active)
% % In normal cases when using L7 segments (except apex and weak)
% % we only need to exclude the size apex segments
% %Seg_active([128, 135, 142, 149, 156, 163]) = [];
% disp('Looping for phase reconstruction...');
% while(true) 
%     iterN = iterN + 1;
%     if mod(iterN,100)==0
%         disp(iterN)
%     end
% 
%     for p = Seg_active  % loop for in-pupil && unlocked segments
%         % determine 6 neighbours surrounding index p (hex-index)        
%         [neighbours_hex, weights_hex] = Find6NeighbrIdx(p, Param, Convert_hex2jk, Convert_jk2hex);
% 
%         % calculate slopes
%         slope_p(1) =  slopes(p, 1);  % rad (pixel / pixel)
%         slope_p(2) =  slopes(p, 2);
%         slopes_hex(:, 1) = slopes(neighbours_hex, 1);
%         slopes_hex(:, 2) = slopes(neighbours_hex, 2);
% 
%         % update wavefront
%         sum_weights = sum(weights_hex);
%         term0 = sum(weights_hex .* W(neighbours_hex)) / sum_weights;
%         term1 = DM_pitch / 4 * ( slope_p(2) * (2*weights_hex(1) - 2*weights_hex(4) + weights_hex(2) ...
%                                                - weights_hex(3) - weights_hex(5) + weights_hex(6))  ... 
%                                + slope_p(1) * sqrt(3) * (weights_hex(2) + weights_hex(3) - weights_hex(5) - weights_hex(6))  );
%         term2 = DM_pitch / 4 * ( (2*weights_hex(1)*slopes_hex(1, 2) - 2*weights_hex(4)*slopes_hex(4, 2) ...
%                                    + weights_hex(2)*slopes_hex(2, 2) - weights_hex(3)*slopes_hex(3, 2) ...
%                                    - weights_hex(5)*slopes_hex(5, 2) + weights_hex(6)*slopes_hex(6, 2) ) ...
%                                   + sqrt(3)*(weights_hex(2)*slopes_hex(2, 1) + weights_hex(3)*slopes_hex(3, 1) ...
%                                    -weights_hex(5)*slopes_hex(5, 1) - weights_hex(6)*slopes_hex(6, 1)) );
%         W_update(p) = term0 + (term1+term2) / sum_weights;         % unit: um
%     end    
%     % end condition - convergency, delta_wf < tolerance
%     error = norm(W_update - W) / norm(W_update+eps);
%     delta_W = [delta_W; error];
%     W = W_update;
%     if error < tol
%         break;
%     end
% end
% % remove phase offset
% W_update = W_update - mean(W_update(:));
% 
% 
% % plot error curve
% fprintf('Finish phase reconstruction, # of iterations is %d. \n', iterN);
% %% set piston, tip, tilt value to IrisAOPTT
% IrisAO_PTT_res = zeros(Param.SegN, 3);
% IrisAO_PTT_res(:, 1) = -W_update;    % z: um
% % IrisAO_PTT_res(:, 1) = 0;   % z: um
% 
% IrisAO_PTT_res(:, 2) = slopes(:, 2) * 1000;  % xgrad, ygrad: mrad
% IrisAO_PTT_res(:, 3) = slopes(:, 1) * 1000;  % xgrad, ygrad: mrad
% % figure(12311); plot(IrisAO_PTT_res(:, 1));
% 
% % add PTT-res to previous IrisAO-PTT
% IrisAO_PTT_res = -IrisAO_PTT_res / 2;  % reflective DM doubles the optical path
% % IrisAO_PTT = IrisAO_PTT_res; 
% IrisAO_PTT = IrisAO_PTT0_1 + IrisAO_PTT_res; 
% IrisAO_PTT([Param.UnreachableSeg],:)=0;
% IrisAO_PTT([Param.LockedSeg],:)=0;
% 
% %% ------------------ post-processing ------------------------------
% 
% 
% % Display and save results
% % SH_name = strrep(SH_name,'SH_IfM\', '');
% % fn = strrep(SH_name, '_', '-');
% % 
% % % display SH image
% % figure(3), set(gcf,'Position',[150 700 400 300]);
% % imagesc(SH_im); axis image; title(['Spot field: ' fn]); colormap hot; colorbar;
% 
% %% calculate & display wavefront - round 1
% wavefront1 = -IrisAO_PTT_res(:, 1)* 2;    % residual aberration     
% wavefront2 = -(IrisAO_PTT(:, 1)) * 2; % collective AO correction
% wavefront_residual = -(IrisAO_PTT(:, 1) - IrisAO_PTT_GT(:, 1)) * 2;
% 
% res_slopes = IrisAO_PTT_res(:, 3:-1:2) * 2 /1000 * SHcam_pitch;
% collective_slopes = IrisAO_PTT(:, 3:-1:2) * 2 /1000 * SHcam_pitch;
% residual_slopes = (IrisAO_PTT(:, 3:-1:2) - IrisAO_PTT_GT(:, 3:-1:2)) * 2 /1000 * SHcam_pitch;
% 
% % SH_im = zeros(1312);
% wavefront1_d = weight_display(SH_im, PosRef, wavefront1, 77, Param.UnreachableSeg, res_slopes); %-slopes/2*SHcam_pitch
% wavefront2_d = weight_display(SH_im, PosRef, wavefront2, 77, Param.UnreachableSeg, collective_slopes);
% wavefrontResidual_d = weight_display(SH_im, PosRef, wavefront_residual, 77, Param.UnreachableSeg, residual_slopes);
% 
% figure(4); set(gcf,'Position',[10 450 1100 300]);
% subplot(1,4,1); imagesc( wavefront1_d - min(wavefront1_d(:)));axis image; colorbar;title(['Measured Aberration: ' num2str(1)]); %colormap jet; 
% subplot(1,4,2);imagesc( wavefront2_d - min(min(wavefront2_d)));axis image; colorbar;title(['Corrective Wavefront: ' num2str(1)]); %colormap jet; 
% subplot(1,4,3);imagesc( wavefrontResidual_d - min(min(wavefrontResidual_d)));axis image; colorbar;title(['Residual Wavefront: ' num2str(1)]); %colormap jet; 
% 
% % wavefrontResidual_d2 = wavefront2_d -wavefront1_d;
% % figure; imagesc( wavefrontResidual_d2 - min(min(wavefrontResidual_d2)));axis image; colorbar;title(['Residual Wavefront: ' num2str(1)]);
% size(wavefront1_d)
% %%
% pupilSize = round((mean(PosRef(158:161, 1)) - mean(PosRef(137:140, 1))) / Param.level * (Param.level + 1/6));
% % centerPos = PosRef(1,:);
% centerPos(1) = round((mean(PosRef(117:121, 1)) + mean(PosRef([99, 101:103], 1))) / 2);
% centerPos(2) = round((mean(PosRef([146, 109, 111, 152], 2)) + mean(PosRef([131, 93, 127, 167], 2))) / 2);
% % set PosRef of Seg1 as the center of the pupil
% % centerPos = [size(wavefront1_d, 1) - PosRef(1, 2), size(wavefront1_d, 2) - PosRef(1, 1)]; 
% figure(4), hold on;
% ang = 0:0.01:2*pi; 
% xp = pupilSize / 2 *cos(ang);
% yp = pupilSize / 2 *sin(ang);
% plot(centerPos(1)+xp, centerPos(2)+yp, 'r');
% hold off;
% centerPos = [size(wavefront1_d, 1) - centerPos(2), size(wavefront1_d, 2) - centerPos(1)]; 
% zernikeCoeff1 = ZernikeDecomposition(rot90(fliplr(wavefront1_d), 3), centerPos, pupilSize);
% zernikeCoeff2 = ZernikeDecomposition(rot90(fliplr(wavefront2_d), 3), centerPos, pupilSize);
% zernikeCoeff3 = ZernikeDecomposition(rot90(fliplr(wavefrontResidual_d), 3), centerPos, pupilSize);
% 
% % zernikeCoeff1 = PTT2ZernikeCoeff(-IrisAO_PTT_res * 2);  % ----------This is a bad mapping
% % zernikeCoeff2 = PTT2ZernikeCoeff( (IrisAO_PTT - IrisAO_PR) * 2);
% figure(5), set(gcf, 'Position', [10 100 800 300]);
% subplot(1,3,1); bar(zernikeCoeff1); title(['Measured Zernike: ' num2str(1)]);
% xlabel('Zernike Mode'); ylabel('Zernike Coefficient (um rms)');
% subplot(1,3,2); bar(zernikeCoeff2); title(['Corrective Zernike: '  num2str(1)]);
% xlabel('Zernike Mode'); ylabel('Zernike Coefficient (um rms)');
% subplot(1,3,3); bar(zernikeCoeff3); title(['Residual Zernike: '  num2str(1)]);
% xlabel('Zernike Mode'); ylabel('Zernike Coefficient (um rms)');
% residual_wavefront1 = wavefrontResidual_d - min(min(wavefrontResidual_d));
% residual_wavefront1_zero = residual_wavefront1 - residual_wavefront1(1,1);
% % figure(1231); imagesc(residual_wavefront1_zero);axis image; colormap jet; hold on;
% mask = poly2mask(centerPos(1)+xp, centerPos(2)+yp, size(residual_wavefront1_zero,1), size(residual_wavefront1_zero,2));
% % grayImage(~mask) = 0; % Blacken outside the mask.
% % imagesc(mask); axis image;
% rms_error1 = std(residual_wavefront1_zero(mask));
% rms_error1_lambda =rms_error1/ 0.488
% 
% %% Remove piston, tip, tilt & defocus
% Param.LockedSeg = [5, 70,  80, 100, 155]; % 20190924: added #80 
% 
% Param.UnreachableSeg = [128, 135, 142, 149, 156, 163]; 
% Param.WeakSeg = [169, 129, 134, 136, 141, 143, 148, 150, 155, 157, 162, 164];
% Param.isRmvPTTD =1;
% %% phase wrapping - round 1
% phase_thres = 0.515; % +- lambda [um]
% phase_wrap = IrisAO_PTT(:, 1);
% while (find(phase_wrap > phase_thres))
%     phase_wrap(phase_wrap > phase_thres) = phase_wrap(phase_wrap > phase_thres) - phase_thres;
% end
% while (find(phase_wrap <- phase_thres))
%     phase_wrap(phase_wrap < -phase_thres) = phase_wrap(phase_wrap < -phase_thres) + phase_thres;
% end
% % IrisAO_PTT(:, 1) = phase_wrap;
% wavefront3 = -(phase_wrap - IrisAO_PTT_GT(:, 1)) * 2; % collective AO correction
% wavefront3_d = weight_display(SH_im, PosRef, wavefront3, Period, Param.UnreachableSeg, collective_slopes);
% figure(4), subplot(1,4,4);
% imagesc( wavefront3_d - min(wavefront3_d(:)));axis image; colorbar;title(['Corrective Wavefront: '  num2str(1)]); %colormap jet; 
% 
% 
% if Param.isRmvPTTD
% 
%     %% update PTT
%     load('ZernikeCoeff2PTT2.mat', 'PTT_mode');
% %     load('ZernikeCoeff2PTT_original_code.mat', 'PTT_mode');
% %     load('ZernikeCoeff2PTT_eye.mat', 'PTT_mode');
%     IrisAO_PTT_res_2 = IrisAO_PTT_res;
%     % IrisAO_PTT_2 =  IrisAO_PTT;
%     IrisAO_PTT_res_2_residual = IrisAO_PTT - IrisAO_PTT_GT;
% 
%     for i=[1:3]
%         IrisAO_PTT_res_2 = IrisAO_PTT_res_2 + squeeze(zernikeCoeff1(i) /2 * PTT_mode(i, :, :));
%         IrisAO_PTT_res_2_residual = IrisAO_PTT_res_2_residual + squeeze(zernikeCoeff3(i) /2 * PTT_mode(i, :, :));
%     end
%     IrisAO_PTT_2 =  IrisAO_PTT0_2 + IrisAO_PTT_res_2;
% 
% 
% 
%     %% for better display
%     IrisAO_PTT_res_2 = CalcLockedSegShift(IrisAO_PTT_res_2, Param.LockedSeg, Param, Convert_hex2jk, Convert_jk2hex);
%     IrisAO_PTT_res_2 = CalcLockedSegShift(IrisAO_PTT_res_2, Param.WeakSeg, Param, Convert_hex2jk, Convert_jk2hex);
%     IrisAO_PTT_2 = IrisAO_PTT0_2 + IrisAO_PTT_res_2;
% 
% 
%     %% calculate & display wavefront - round 2
%     wavefront1 = -IrisAO_PTT_res_2(:, 1) * 2;         
%     wavefront2 = -(IrisAO_PTT_2(:, 1)) * 2; % collective AO correction
%     wavefront_residual = -(IrisAO_PTT_2(:, 1) - IrisAO_PTT_GT(:, 1)) * 2;
%     wavefront_gt = -IrisAO_PTT_GT(:, 1) * 2;
% 
%     res_slopes = IrisAO_PTT_res_2(:, 3:-1:2) * 2 /1000 * SHcam_pitch;
%     collective_slopes = (IrisAO_PTT_2(:, 3:-1:2)) * 2 /1000 * SHcam_pitch;
%     residual_slopes = (IrisAO_PTT_2(:, 3:-1:2) - IrisAO_PTT_GT(:, 3:-1:2)) * 2 /1000 * SHcam_pitch;
%     gt_slopes = IrisAO_PTT_GT(:, 3:-1:2) * 2 /1000 * SHcam_pitch;
% 
%     wavefront1_d = weight_display(SH_im, PosRef, wavefront1, 77, Param.UnreachableSeg, res_slopes); %-slopes/2*SHcam_pitch
%     wavefront2_d = weight_display(SH_im, PosRef, wavefront2, 77, Param.UnreachableSeg, collective_slopes);
%     wavefrontResidual_d = weight_display(SH_im, PosRef, wavefront_residual, 77, Param.UnreachableSeg, residual_slopes);
%     wavefrontGT_d = weight_display(SH_im, PosRef, wavefront_gt, 77, Param.UnreachableSeg, gt_slopes);
%     residual_wavefront2 = wavefrontResidual_d - min(min(wavefrontResidual_d));
%     residual_wavefront2_zero = residual_wavefront2 - residual_wavefront2(1,1);
%     % figure(1231); imagesc(residual_wavefront1_zero);axis image; colormap jet; hold on;
%     mask = poly2mask(centerPos(1)+xp, centerPos(2)+yp, size(residual_wavefront2_zero,1), size(residual_wavefront2_zero,2));
%     % grayImage(~mask) = 0; % Blacken outside the mask.
%     % imagesc(mask); axis image;
%     rms_error2 = std(residual_wavefront2_zero(mask));
%     rms_error2_lambda = rms_error2/ 0.488
% 
%     figure(6); set(gcf,'Position',[810 450 1100 300]); 
% %     [xx, yy] = meshgrid(1:size(wavefront1_d, 1), 1:size(wavefront1_d, 2));
% %     subplot(1,3,1); surf(xx, yy, wavefront1_d,'EdgeColor','none');
% %     subplot(1,3,2); surf(xx, yy, wavefront2_d,'EdgeColor','none');
%     subplot(1,5,1); imagesc(wavefront1_d - min(wavefront1_d(:)));axis image; xticks([]); yticks([]); colorbar;title(['Measured Aberration: ' num2str(1)]); colormap jet; 
%     subplot(1,5,2);imagesc(wavefront2_d - min(wavefront2_d(:)));axis image; colorbar;title(['Corrective Wavefront: ' num2str(1)]); colormap jet; 
%     subplot(1,5,3);imagesc(wavefrontResidual_d - min(wavefrontResidual_d(:)));axis image; colorbar;title(['Residual Wavefront: ' num2str(rms_error2_lambda)]); colormap jet; 
% 
%     preserveOTF = false;
%     if (preserveOTF)
%         alpha_data = ones(size(wavefront2_d));
%         for y = 1:size(wavefront2_d, 1)
%             for x =1:size(wavefront2_d, 2)
%                 r2 = (y - centerPos(2))^2 + (x -centerPos(1))^2;
%                 if r2 > (pupilSize/2)^2
%                     alpha_data(y,x) = 0;
%                     wavefront2_d(y,x) = wavefront2_d(1,1);
%                 elseif r2 > (pupilSize/2 - 5)^2
%                     alpha_data(y,x) = cos((r2-(pupilSize/2-5)^2)/((pupilSize/2)^2-(pupilSize/2-5)^2) * pi/2)^2;
%                 end
%             end
%         end
%         figure(8), h = imshow(-wavefront2_d+min(wavefront2_d(:)),[]);colormap(parula);colorbar;
%         set(h, 'AlphaData', alpha_data);
%         pause();
%         return;
%     end
% 
%     %% zernike coeff: remove mode #1-3,5 & display
%     zernikeCoeff1([1:3, 5]) = 0; % remove 1-3, 5
%     zernikeCoeff2([1:3, 5]) = 0; % remove 1-3, 5
%     zernikeCoeff3([1:3, 5]) = 0; % remove 1-3, 5
% 
%     figure(7), set(gcf, 'Position', [810 100 800 300]);
%     subplot(1,3,1); bar(zernikeCoeff1); title(['Measured Zernike: ' num2str(1)]);
%     xlabel('Zernike Mode'); ylabel('Zernike Coefficient (um rms)');
%     subplot(1,3,2); bar(zernikeCoeff2); title(['Corrective Zernike: ' num2str(1)]);
%     xlabel('Zernike Mode'); ylabel('Zernike Coefficient (um rms)');
%     subplot(1,3,3); bar(zernikeCoeff3); title(['Residual Zernike: ' num2str(1)]);
%     xlabel('Zernike Mode'); ylabel('Zernike Coefficient (um rms)');
% 
%     %% phase wrapping - round 2
%     phase_wrap_2 = IrisAO_PTT_2(:, 1);
%     while (find(phase_wrap_2 > phase_thres))
%         phase_wrap_2(phase_wrap_2 > phase_thres) = phase_wrap_2(phase_wrap_2 > phase_thres) - phase_thres;
%     end
%     while (find(phase_wrap_2 <- phase_thres))
%         phase_wrap_2(phase_wrap_2 < -phase_thres) = phase_wrap_2(phase_wrap_2 < -phase_thres) + phase_thres;
%     end
%     disp(123)
%     % res_slopes
%     % wavefront3
% %     IrisAO_PTT_2(:, 1) = phase_wrap_2;
%     wavefront3 = -(phase_wrap_2 - IrisAO_PTT_GT(:, 1)) * 2; % collective AO correction
%     wavefront3_d = weight_display(SH_im, PosRef, wavefront3, 77, Param.UnreachableSeg, res_slopes);
%     % wavefront3_d = weight_display(SH_im, PosRef, wavefront3, 77, [], res_slopes);
% 
%     figure(6), subplot(1,5,4);
% %     s = surf(xx, yy, wavefront3_d,'EdgeColor','none');
%     imagesc( wavefront3_d - min(wavefront3_d(:)));axis image; colorbar;title(['Corrective Wavefront: ' 1]);  %colormap jet; 
%     figure(6), subplot(1,5,5);
% %     s = surf(xx, yy, wavefront3_d,'EdgeColor','none');
%     imagesc( wavefrontGT_d - min(wavefrontGT_d(:)));axis image; colorbar;title(['gt Wavefront: ' 1]); %caxis([0.4,1]) %colormap jet;  
% 
% end
% %% update PTT - round 3 
% IrisAO_PTT_wrap = IrisAO_PTT;
% IrisAO_PTT_wrap(:, 1) = phase_wrap;
% if Param.isRmvPTTD
%     IrisAO_PTT_2_wrap = IrisAO_PTT_2;
%     IrisAO_PTT_2_wrap(:, 1) = phase_wrap_2;
% end
% %% set all outlier's PTT = 0
% % IrisAO_PTT(Param.Outlier_hex,:)=0;
% IrisAO_PTT_DM_wrap = IrisAO_PTT_wrap;
% IrisAO_PTT_DM_wrap(Param.Outlier_hex,:)=0;
% % IrisAO_PTT_2(Param.Outlier_hex,:)=0;
% IrisAO_PTT_2_DM_wrap = IrisAO_PTT_2_wrap;
% IrisAO_PTT_2_DM_wrap(Param.Outlier_hex,:)=0;
% %%
% 
% function [shift] = single_shift_finder(SH_im, pos_init, PosRef, Period, isDisplay, Sig, ItN)
% % 
% halfSize = floor(Period/2);
% % if halfSize > size(SH_im,2)-pos_init(1) ||halfSize > size(SH_im,1)-pos_init(2) || halfSize > pos_init(2)-1 || halfSize > pos_init(1)-1
% % 	shift = [nan, nan];
% %     return
% % end
% 
% rect = SH_im(pos_init(2)-halfSize:pos_init(2)+halfSize,pos_init(1)-halfSize:pos_init(1)+halfSize);
% % if pos_init(1)<200 && pos_init(2)<200 
% %     shift = [nan, nan];
% % else
% 
%     [x0,y0] = peak_finder(rect);
%     % gausss fit
% 	% disp([x0,y0]);
% 
%     ChosenPixels = 2;
%     if pos_init(2) + y0-ChosenPixels < 1 ||pos_init(2)+y0+ChosenPixels > size(SH_im,1) || pos_init(1) + x0-ChosenPixels <1 || pos_init(1)+x0+ChosenPixels > size(SH_im,1)
%         % shift = [pos_init(1)+x0, pos_init(2)+y0];
% 	    shift = [x0, y0];
%         % shift = [nan, nan];
%         return
%     end
%     rect_new = SH_im(pos_init(2) + y0-ChosenPixels:pos_init(2)+y0+ChosenPixels, pos_init(1) + x0-ChosenPixels:pos_init(1)+x0+ChosenPixels);
%     [x1, y1] = gaussian_mask_fit(rect_new, Sig, ItN);
%     shift = [0,0];
%     shift(1) = pos_init(1) +x0+x1-PosRef(1);
%     shift(2) = pos_init(2) +y0+y1-PosRef(2);
% % end
% 
% 
% 
% if isDisplay
%     figure(101);hold on; 
%     %scatter(PosRef(1), PosRef(2), 'r*');         % original PosRef
%     if ~isnan(shift(1))
%         rectangle('Position',[pos_init(1)-halfSize pos_init(2)-halfSize halfSize*2 halfSize*2]);
%         scatter(PosRef(1)+shift(1), PosRef(2)+shift(2), 'g'); % Shifted Peakpoint     
%     end
% 
% end
% end
% 
% function [x0,y0] = peak_finder(Data)
% N = size(Data,1);
% [~, id] = max(Data(:));
% x0 = ceil(id/N);
% y0 = id-N*(x0-1);
% x0 = x0-(N+1)/2;
% y0 = y0-(N+1)/2;
% end
% 
% function [x,y] = gaussian_mask_fit(RawData, Sig, ItN)
% % Sig = Param.Sig;
% % ItN = Param.ItN;
% 
% x = 0;
% y = 0;
% N = size(RawData,1);
% [X,Y] = meshgrid(-(N-1)/2:1:(N-1)/2,-(N-1)/2:1:(N-1)/2);
% for it = 1:ItN
%     PSF = exp(-1/2/Sig^2*((X-x).^2+(Y-y).^2));
%     x = sum(sum(X.*RawData.*PSF)) / (sum(sum(RawData.*PSF))+eps);
%     y = sum(sum(Y.*RawData.*PSF)) / (sum(sum(RawData.*PSF))+eps);
% end
% 
% end

% Define the 2D Gaussian function
function F = gaussian2D(params, X)
    % Unpack parameters
    A = params(1); % Amplitude
    x0 = params(2); % Center in X
    y0 = params(3); % Center in Y
    sigma_x = params(4); % Std dev in X
    sigma_y = params(5); % Std dev in Y
    % C = params(6); % Offset

    % Create a 2D Gaussian function
    x = X(:, 1);
    y = X(:, 2);
    
    % F = A * exp(-((x - x0).^2 / (2 * sigma_x^2) + (y - y0).^2 / (2 * sigma_y^2))) + C;
    F = A * exp(-((x - x0).^2 / (2 * sigma_x^2) + (y - y0).^2 / (2 * sigma_y^2))) + 0;

end
