% =============================================================
%    Fourier-plane mask demo   (±1 orders of three gratings)
%    — all scaling driven by NA, λ, pixel size, magnification —
% =============================================================
close all; clear; clc;

%% ------------------  Parameter Block  ------------------------
N              = 512;        % simulation grid (square)
lambda         = 0.515e-6;   % wavelength in vacuum  [m]
NA             = 1.1;        % numerical aperture
pixelSize_cam  = 6.5e-6;     % sensor pixel pitch    [m]
M              = 75.6303;    % lateral magnification  (sample→camera)

period_slm =3.97;
% period_slm =5.95;

period_slm_um = period_slm *8.2;
mag_Slm_Camera = (125/150)*(300/175)*(75/85);
% rho = ((period_slm*mag_Slm_Camera)^-1)./((size(Image_temp,1)*6.5)^-1)*ones(4,1);
period = period_slm_um*mag_Slm_Camera / 6.5;
% period         = 6;          % grating period [pixels]  (simulation coord.)
maskrad        = 3;          % radius of each ±1 order mask [pixels]
angles_list    = -[94.82, 94.82-120, 94.82-240];  % grating orientations

%% ------------------  Derived constants -----------------------
% frequency increment on the camera (cycles/m per FFT pixel)
deltaF_cam = 1 / (N * pixelSize_cam);
% res_rayleigh = 0.61 * lambda / NA;          % [m]

% NA cut-off in cycles/m (object) and camera
% fc_obj = NA / (0.61 * lambda);    % Rayleigh's limit - object plane
fc_obj = 2*NA / (1 * lambda);       % Abbe's limit - object plane

fc_cam = fc_obj / M;         % camera plane  (what the FFT sees)

% radius of NA circle in pixels
kCutPix = round(fc_cam / deltaF_cam);

fprintf('NA cut-off (object) : %.3f cycles/µm\n',  fc_obj*1e-6);
fprintf('NA cut-off (camera) : %.3f cycles/µm\n',  fc_cam*1e-6);
fprintf('NA radius in FFT px : %d px\n\n',         kCutPix);

%% ------------------  Build spatial grid ----------------------
[x, y] = meshgrid( (0:N-1) - N/2, (0:N-1) - N/2 );

%% ------------------  Composite binary-phase field ------------
E_tot = ones(N);
for jj = 1:3
    th     = angles_list(jj);
    coords = x*cosd(th) + y*sind(th);
    base_phase = pi * (mod(coords,period) < period/2);   % 0–π stripes
    E_tot = E_tot + exp(1i*base_phase);
end

figure(1);
subplot(1,2,1); imagesc(abs(E_tot));  axis image; colorbar; colormap jet;
title('|E_{tot}|');
subplot(1,2,2); imagesc(angle(E_tot)); axis image; colorbar; colormap jet;
title('arg(E_{tot})');

%% ------------------  Fourier transform -----------------------
Spec = fftshift( fft2(E_tot) );

%% ------------------  Mask: keep ±1 orders --------------------
mask = false(N);
center = N/2 + 1;

% ±1 orders for each orientation
f0 = N/period;   % FFT-bin index of first order
[Xg, Yg] = meshgrid(1:N, 1:N);

for th = angles_list
    u = round(f0 * cosd(th));
    v = round(f0 * sind(th));
    for sgn = [+1, -1]
        rr = center + sgn*v;
        cc = center + sgn*u;
        mask = mask | ((Yg-rr).^2 + (Xg-cc).^2 <= maskrad^2);
    end
end

% ------------------  Add circular NA pupil -------------------
NAMask  = ((Xg-center).^2 + (Yg-center).^2) <= kCutPix^2;
% mask    = mask & NAMask;   % keep only ±1 orders *inside* NA circle
mask    =  NAMask;   % keep only ±1 orders *inside* NA circle

figure(2); imagesc(mask); axis image;
title('Composite mask:  ±1 orders  ∩  NA pupil');
ft_profile = abs(Spec);
ft_profile = ft_profile ./ max(ft_profile(:));
%%
figure(3); imagesc(ft_profile); axis image off; %caxis([3e-2, 1e-1]);
set(gca,'ColorScale','log'); 
colormap gray;colorbar; hold on;
caxis([3e-3, 1e-1]);
viscircles([center,center], kCutPix, 'LineStyle','--','Color','w');
title('|FT| with NA circle (camera plane)');

%% ------------------  Filter, inverse FT ----------------------
Spec_filt = Spec .* mask;

% If you want to flip the phase of certain gratings:
to_flip = [];   % e.g. [1 2] to flip the first two gratings
for k = to_flip
    thk = angles_list(k);
    uk  = round(f0 * cosd(thk));
    vk  = round(f0 * sind(thk));
    for sgn = [+1, -1]
        rr = center + sgn*vk;
        cc = center + sgn*uk;
        Spec_filt(rr,cc) = -Spec_filt(rr,cc);   % multiply by –1
    end
end

E_img = ifft2( ifftshift(Spec_filt) );
I_img = abs(E_img).^2;

%% ------------------  Display focal-plane intensity -----------
figure(7); imagesc(I_img); axis image off; colormap gray; colorbar;
title('Focal‐plane intensity after ±1 orders of three gratings');
%% ----- Compute ±1 order radius ratio vs NA radius -----
% Continuous (ideal) first-order radius in FFT pixels
r1_cont = f0;                       % since radius = N/period in pixels
ratio_cont = r1_cont / kCutPix;     % continuous ratio

% Discrete (rounded) radius from mask centers for each grating angle
r1_disc = zeros(numel(angles_list), 1);
cnt = 0;
for th = angles_list
    cnt = cnt + 1;
    u = round(f0 * cosd(th));
    v = round(f0 * sind(th));
    r1_disc(cnt) = hypot(u, v);     % sqrt(u^2 + v^2)
end
ratio_disc_all = r1_disc / kCutPix;
ratio_disc_mean = mean(ratio_disc_all);

% Print to console
fprintf('--- ±1 Order Radius vs NA ---\n');
fprintf('NA radius (kCutPix)           : %d px\n', kCutPix);
fprintf('First-order radius (continuous): %.3f px\n', r1_cont);
fprintf('Radius ratio (continuous)      : %.6f\n', ratio_cont);
for ii = 1:numel(angles_list)
    fprintf('Angle %7.2f°: r_disc = %4d px  |  ratio_disc = %.6f\n', ...
        angles_list(ii), round(r1_disc(ii)), ratio_disc_all(ii));
end
fprintf('Mean discrete ratio across angles: %.6f\n\n', ratio_disc_mean);
