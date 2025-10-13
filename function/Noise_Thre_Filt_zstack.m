function [out_filt_zstack, max_ind, out_max, noise_bg] = Noise_Thre_Filt_zstack(raw_img_stack, thr_para, kernel_size, pix_cnt)
    z_length = size(raw_img_stack,3);

    % find mean, std intensity
    filt_img_zstack = raw_img_stack;
    max_filt_zstack = zeros(1, z_length);
    for i=1:z_length
        if kernel_size>0
            filt_img_zstack(:,:,i) = medfilt2(raw_img_stack(:,:,i), [kernel_size, kernel_size], 'symmetric');
        else
            filt_img_zstack(:,:,i) = raw_img_stack(:,:,i);
        end
        max_filt_zstack(i) = max(max(filt_img_zstack(:,:,i)));
    end
    [~, max_ind]= max(max_filt_zstack);
    filt_img_peak = filt_img_zstack(:,:,max_ind);
    mean_intensity_peak = mean(filt_img_peak(:));
    std_inetnsity_peak = std(filt_img_peak(:));

    threshold_value = mean_intensity_peak + thr_para.* std_inetnsity_peak; % use local threshold 
    Binary_Im = filt_img_peak < threshold_value;  %binary image under threshold
    figure(111)
    imagesc(Binary_Im); axis image; colormap hot; title([max_ind]);
    Noise_avg = mean(mean(filt_img_peak(Binary_Im)))
    
    out_filt_zstack = filt_img_zstack - Noise_avg;
    out_filt_mean_zstack = zeros(1, z_length); 
    out_filt_max_zstack = zeros(1, z_length);

    [M, img_ind] = max(filt_img_peak(:));
    [max_x, max_y] = ind2sub(size(filt_img_peak,2), img_ind);
    % disp([M, img_ind])
    % figure; imagesc(filt_img_peak);
    for i=1:z_length
        result_zstack = out_filt_zstack(:,:,i);
%          [M, img_ind] = max(result_zstack(:));
%         [max_x, max_y] = ind2sub(size(result_zstack,2), img_ind);
        binary_result_filt = result_zstack>threshold_value;
        result_zstack_temp = result_zstack.*binary_result_filt;
%         imagesc(Binary_Im); axis image; colormap hot;
        out_filt_mean_zstack(i) = mean(mean(result_zstack(max_x-pix_cnt:max_x+pix_cnt,max_y-pix_cnt:max_y+pix_cnt))); 
        out_filt_max_zstack(i) = max(max(result_zstack(:)));
    end

    out_max = out_filt_mean_zstack;
%     out_max = out_filt_max_zstack;
    noise_bg = Noise_avg;
end

