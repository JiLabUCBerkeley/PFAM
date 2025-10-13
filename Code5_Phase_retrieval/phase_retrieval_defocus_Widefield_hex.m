
function [Phase, OpticalPath, PupilEst, PupilMask] = ...
    phase_retrieval_defocus_Widefield_hex(ImStack_name, IrisAORef_name, defocus)
%% =======================================
% load in "syscor-new.tif", "syscor-old.txt"
% save to "syscor-new.txt", "syscor-new.mat"

if nargin < 3
    defocus = 0;
end

Center =[66, 60];            % (y, x) coordinates of the bead center; default [66, 60]
ZernikeModeN = 300;
Param.Lamd=0.515; % default green fluorescence 515
Param.NA= 1.1; %1.1, 0.45, 0.25
Param.RI=1.33; %1.33, 1
Param.Pixel=0.086; % default WF 0.086 (1.1), 0.225, 0.35
Param.ImageSize=61;         % cropped region with centered bead default 61
Param.ImageEdge=20;
Param.ItN=45; 
Param.Unwrap=1;
Param.FixAmp=1;
Param.BackgroundCoef=0.998;
Param.BackgroundOffset=200;  % background
Param.IntegrationN=1;
Param.RemoveTilt=1;

load([IrisAORef_name '.mat'], 'IrisAO_Modal');
IrisAO_Modal = reshape(IrisAO_Modal, numel(IrisAO_Modal), 1); % force to be Nx1-dimension
IrisAORef = IrisAO_Modal;

ImStack = img_get_experimental(Param,ImStack_name,Center);
disp((length(ImStack(1,1,:))-1))
% 25X 1.1NA [0.7~-0.7] @0.05 - VSIM input, real value needs to time 4
% partial 21 stack
Param.ZDefocus = (-2:4/(length(ImStack(1,1,:))-1):2)*(-1)+defocus; 

% 10X 0.45NA [2.4~-2.4] @0.2
% Param.ZDefocus = (-8:16/(length(ImStack(1,1,:))-1):8)*(-1)+defocus;    

% 5X 0.25NA [2.4~-2.4] @0.2
% Param.ZDefocus = (-40:80/(length(ImStack(1,1,:))-1):40)*(-1)+defocus;    

[PupilEst, PupilMask] = phase_retrieval(Param,ImStack);

Phase = QualityGuidedUnwrap2D(PupilEst, PupilMask,...
    [(size(PupilEst,1)+1)/2 (size(PupilEst,1)+1)/2]);

%% phase to Zernikecoeff
OpticalPath = Phase / 2 / pi * Param.Lamd;
wvf_rms = norm(OpticalPath) / numel(OpticalPath);
figure(100); imagesc(OpticalPath); axis image; colorbar;
title(sprintf('Wavefront error [um], rms = %.2e', wvf_rms));
IrisAO_Modal = Phase2IrisAO(OpticalPath, ZernikeModeN); % 21 rms-value
IrisAO_Modal = reshape(IrisAO_Modal, numel(IrisAO_Modal), 1); % force to be Nx1-dimension

%% display ZernikeCoeff
IrisAO_Modal(1:3) = IrisAO_Modal(1:3)*0; 
figure(3); bar(IrisAO_Modal(1:55)); title('Zernike Coeff');
xlabel('Modes'); ylabel('Amplitude RMS [um]');
IrisAO_Modal(5) = IrisAO_Modal(5)*0; 

%% save value
IrisAO_Modal = IrisAO_Modal / 2 ; % DM doubles the optical path
% IrisAO_Modal = IrisAO_Modal / 4 ; % DM doubles the optical path 201907223 test /2
IrisAO_Modal = -IrisAO_Modal + IrisAORef;
save(['data/' ImStack_name '.mat'],'IrisAO_Modal');
figure(4); bar(IrisAO_Modal(1:55)); title('Total Zernike Coeff');


%% convert 1-55 modal coeffs to PTT
orderN = 55;
IrisAO_PTT_LOrder = zeros(169, 3);
IrisAO_PTT_HOrder = zeros(169, 3);
load('data/ZernikeCoeff2PTT2.mat', 'PTT_mode');
for i = 1:21
    IrisAO_PTT_LOrder = IrisAO_PTT_LOrder + squeeze(PTT_mode(i,:,:) * IrisAO_Modal(i));
end
for i = 22:orderN
    IrisAO_PTT_HOrder = IrisAO_PTT_HOrder + squeeze(PTT_mode(i,:,:) * IrisAO_Modal(i));
end
IrisAO_PTT = IrisAO_PTT_LOrder + IrisAO_PTT_HOrder;

save(['data/' ImStack_name '.mat'], 'IrisAO_PTT', 'IrisAO_PTT_LOrder', 'IrisAO_PTT_HOrder', '-append');
clear PTT_mode;

end



function ImStack=img_get_experimental(Param,FileName, Center)
img = imstackread([FileName '.tif']); % function not found
img = reshape(img,[size(img,1) size(img,2),...
    Param.IntegrationN size(img,3)/Param.IntegrationN]); % y, x, I, z/I
img = squeeze(sum(img,3)); % Integration of successive frames: y, x, z/I
img = double(img); 

% find center point with highest intensity
Img0 = img(:,:, (size(img, 3)+1)/2);
% only llok at the central region
[w, h] = size(Img0);
Img0_central = Img0(round(w/4):round(w*3/4), round(h/4):round(h*3/4)); % 20210408, added by YH
%
[mI, mY] = max(Img0_central); % row vector
[~, mx] = max(mI);
my = mY(mx);
Center = [my + round(w/4) - 1, mx + round(h/4) - 1];
disp('The center point coordination is:');
disp(Center);
% find background intensity
threshold = median(Img0(:)) * 2;
Param.BackgroundOffset = mean2(Img0(Img0<threshold));
disp('The background intensity is:');
disp(Param.BackgroundOffset);

y0 = Center(1)-(Param.ImageSize(1)-1)/2; 
y1 = Center(1)+(Param.ImageSize(1)-1)/2;
x0 = Center(2)-(Param.ImageSize(1)-1)/2;
x1 = Center(2)+(Param.ImageSize(1)-1)/2;

for ii=1:size(img,3)
    ImStack(:,:,ii)=img(y0:y1,x0:x1,ii); % crop image
end

for ii=1:size(ImStack,3)
    ImStack(:,:,ii)=rot90(ImStack(:,:,ii), 1); % 90 degree clockwise rotation
end

end


function [PupilEst, PupilMask, residualFieldEst] = phase_retrieval(Param,ImStack)

%% -------------------------------------------------------------------
%% Phase Retrieval: F(k)=|F(k)|exp(j*phi(k))=fft(f(x)), 
%%                  given |F(k)|, find f(x), or equavelently phi(k)
%% G-Saxton Algorithm: initilize f_e(x) 
%%                     loop: 
%%                          F_e(k)=fft(f_e(x)) --> 
%%                          phi_e(k)=phase(F_e(k))
%%                          F_e(k)=|F(k)|*exp(j*phi_e(k)) --> 
%%                          f_e(x)=F^-1(F_e(k))
%%                      end
%% In our case, Field<->f(x), Pupil<->F(k), |Field| is calculated 
%% from image stack, Pupil is initialized to be the OTF. The GS process
%% should change to F(k)-> f(x)-> |f|exp(j*phase(f(x)))-> F(k)-> phase(F(k))
%% -------------------------------------------------------------------

ImStack=ImStack-Param.BackgroundOffset; % background subtraction
ImStack=ImStack.*(ImStack>0);
Field=ImStack.^(1/2); % Image plane: intensity to amplitude
Field=Field-Param.BackgroundCoef*mean(Field(:)); % background subtraction
Field=Field.*(Field>0);
Field=padarray(Field,[Param.ImageEdge,Param.ImageEdge,0]);

ImgSize=size(Field,1);
dk=2*pi/ImgSize/Param.Pixel; % unit spatial frequency ~0.5um^-1
kx=(-(ImgSize-1)/2:1:(ImgSize-1)/2)*dk; % spatial frequency
[kx, ky]=meshgrid(kx,kx);
kz=sqrt((2*pi/Param.Lamd*Param.RI)^2-kx.^2-ky.^2); % ?? kx^2+ky^2+kz^2 = (1.33*2*pi\lambda)^2
PupilMask=(kx.^2+ky.^2)<=(2*pi/Param.Lamd*Param.NA)^2; % Pupil plane: 2NA/lambda, WHY *PI???
PupilN=round(2*pi/Param.Lamd*Param.NA/dk)*2+1; % cutoff frequency in the unit of pixel

PupilEst=PupilMask;
for ii=1:Param.ItN % G-S iteration
    for jj=1:length(Param.ZDefocus)
        FieldEst(:,:,jj)=fftshift(fft2(ifftshift(... 
            PupilEst.*exp(1j*kz*Param.ZDefocus(jj)))));% Fe = ifft(Pe)
    end
    FieldEst=Field.*exp(1j*angle(FieldEst)); % Fe = |F|*exp(j*Phase(Fe))
    PupilEst=0;
    for jj=1:length(Param.ZDefocus)
        PupilEst=PupilEst+fftshift(ifft2(ifftshift(FieldEst(:,:,jj))))...
            .*exp(-1j*kz*Param.ZDefocus(jj)); % Pe = fft(Fe)
    end
    PupilEst=PupilEst/length(Param.ZDefocus);   
    if Param.FixAmp==1
        PupilEst=PupilMask.*exp(1j*angle(PupilEst)); % Mask constraint?
    else
        PupilEst=PupilEst.*PupilMask;
    end
    PupilEstPhase=angle(PupilEst); % phi_e = phase(Pe)
    if Param.RemoveTilt==1
        PupilEstPhase = PupilEstPhase - ...
            sum(PupilEstPhase(:).*kx(:).*PupilMask(:)) / ...
            sum(kx(:).*kx(:).*PupilMask(:)) *kx .*PupilMask;
        PupilEstPhase = PupilEstPhase - ...
            sum(PupilEstPhase(:).*ky(:).*PupilMask(:)) / ...
            sum(ky(:).*ky(:).*PupilMask(:)) *ky .*PupilMask;
    end
%     figure(2);subplot(1,2,1);imagesc(abs(PupilEst));axis image; colorbar; title('pupil mask');
%     figure(2);subplot(1,2,2);imagesc(PupilEstPhase);axis image; colorbar; title('estimated phase');
end

PupilEst = PupilMask.*exp(1j*PupilEstPhase);
PupilMask = PupilMask((ImgSize-PupilN)/2+2:end-(ImgSize-PupilN)/2-1, ...
    (ImgSize-PupilN)/2+2:end-(ImgSize-PupilN)/2-1);
PupilEst = PupilEst((ImgSize-PupilN)/2+2:end-(ImgSize-PupilN)/2-1, ...
    (ImgSize-PupilN)/2+2:end-(ImgSize-PupilN)/2-1);

end


function PhaseZernikeCoef = Phase2IrisAO(Phase, ZernikeModeN)
ImSize = size(Phase,1);
PhaseZernikeCoef = zeros(1, ZernikeModeN);
for ii=1:ZernikeModeN
    [Mode, ~] = zernike_fun(ii,ImSize);                                % Mode: 51x51 wavefront for each Zernike mode
    PhaseZernikeCoef(ii) = sum(Phase(:).*Mode(:)) / sum(Mode(:).^2);   % Zernike_coef*mode=phase --> Zer_coef=phase/mode 
end
end

