function PTT_mode = ZernikeCoeff2PTT(orderN, coeff, pupilSize, Period);

%% Parameters based on imaging system (0.50mm /6.5um = 77)
imSize = 1400;
Period = 77;    % the most important value! MUST be same with the real system!
orderN = 55;
coeff = 0.2;   % to avoid too large wavefront distortion

%% lens-array sampling
N_Level = 7;
SegN = N_Level*(N_Level+1)*3 + 1;
% outliers = [128, 135, 142, 149, 156, 163, 5, 70, 100, 155];
outliers = [128, 135, 142, 149, 156, 163];
% outliers = [];

Region = [imSize/2 - Period*N_Level/2*sqrt(3), imSize/2 - Period*N_Level; 
          imSize/2 + Period*N_Level/2*sqrt(3), imSize/2 + Period*N_Level];
pupilSize = round(Period * (N_Level*sqrt(3) + 1/sqrt(3)));
border = (imSize - pupilSize)/ 2;

step_vec = [-sqrt(3)/2,1/2; 0,1; sqrt(3)/2,1/2;...
            sqrt(3)/2,-1/2; 0,-1; -sqrt(3)/2,-1/2];
        
PosInt = zeros(SegN, 2);
PosInt(1,:) = (Region(1,:)+ Region(2,:)) /2; % level 0
count = 1;
for L = 1:N_Level
    count = count + 1;
    PosInt(count,:) = PosInt(1,:) - [0, Period*L]; % start point for level L
    for i = 1:6 % six orientations
        for j = 1:L % l segments
            if i+j == 6+L
                break;
            end
            count = count + 1;
            PosInt(count,:) = PosInt(count-1,:) + Period * step_vec(i,:); 
        end
    end
end
PosInt = round(PosInt);

   
%% generate zernike mode wavefront
PTT_mode = zeros(orderN, SegN, 3);
for order = 1: orderN
    fprintf('Calculating order #%d...\n', order);
    
    [mode, ~] = zernike_fun(order, pupilSize);
    mode = mode * coeff;
    mode = padarray(mode, [border border]);
    % figure(2), imshow(imresize(mode, 0.2),[])%, colormap parula; colorbar; title('wavefront');
    % mode = mode * (-1); % rms = 1.0, on deformable mirror
    mode = rot90(mode, 1);
    mode = fliplr(mode);

    % figure(2), hold on, scatter(PosInt(:,1)*0.2, PosInt(:,2)*0.2, 'k'), hold off;

    %% calculate gradients
    piston = zeros(SegN, 1);
    gradx = zeros(SegN, 1);
    grady = zeros(SegN, 1);
    rect_Size = round(Period * 0.45);
    active_Seg = 1 : SegN;
    active_Seg(outliers) = [];
    for i = active_Seg
        rect = mode( PosInt(i,2) - rect_Size : PosInt(i,2) + rect_Size, ...
            PosInt(i,1) - rect_Size : PosInt(i,1) + rect_Size );

        [xx, yy] = meshgrid(-rect_Size:rect_Size);       
        % remove out-of pupil points
        xx = xx + PosInt(i, 1) - imSize/2;
        yy = yy + PosInt(i, 2) - imSize/2;
        outofpupil = (xx.^2+yy.^2 >= (pupilSize/2-2)^2);% Qinrong corrected 20200801 from (-1) to (-2)
        outofpupil = outofpupil(:);
        xx = xx - PosInt(i, 1) + imSize/2;
        yy = yy - PosInt(i, 2) + imSize/2;
        xx = xx(:); 
        xx(outofpupil) = [];
        yy = yy(:);
        yy(outofpupil) = [];
        rect = rect(:);
        rect(outofpupil) = [];
        % fit with f=ax+by+c
        fo = fit([xx, yy], rect, 'poly11');
        piston(i) = fo.p00;
        gradx(i) = fo.p10;
        grady(i) = fo.p01;
    %         if sum(outofpupil(:))
    %             fprintf('Seg %d, # of out-of-pupil pixels: %d\n', i, sum(outofpupil(:)));
    %         end
    end

    %% calculate shifts
    lens_pitch = 0.5; % lensarray dimension, milimeter
    gradx = -gradx * Period / lens_pitch; % um/pixel --> um/mm(mrad) 
    grady = -grady * Period / lens_pitch; % um/pixel --> um/mm(mrad)
    shift(order, :) = [gradx; grady] * Period; % 55*338

    %% save PTT values
    PTT_mode(order, :, :) = [piston, grady, gradx];
    
end
disp(outliers);
%% save convertion to .mat
PTT_mode = PTT_mode / coeff;
save('data/ZernikeCoeff2PTT_nooutlier.mat', 'PTT_mode', 'shift');   
   
end