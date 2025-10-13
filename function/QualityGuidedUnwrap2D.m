function im_unwrapped = QualityGuidedUnwrap2D(IM, im_mask, startingpoint)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% QualityGuidedUnwrap2D implements 2D quality guided path following phase
% unwrapping algorithm.
%
% Technique adapted from:
% D. C. Ghiglia and M. D. Pritt, Two-Dimensional Phase Unwrapping:
% Theory, Algorithms and Software. New York: Wiley-Interscience, 1998.
% 
% Inputs: 1. Complex image in .mat double format
%         2. Binary mask (optional)          
% Outputs: 1. Unwrapped phase image
%          2. Phase quality map
%
% This code can easily be extended for 3D phase unwrapping.
% Posted by Bruce Spottiswoode on 22 December 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

im_mag=abs(IM);                             %Magnitude image
im_phase=angle(IM);                         %Phase image
im_unwrapped=zeros(size(IM));               %Zero starting matrix for unwrapped phase
adjoin=zeros(size(IM));                     %Zero starting matrix for adjoin matrix
unwrapped_binary=zeros(size(IM));           %Binary image to mark unwrapped pixels

%% Calculate phase quality map
im_phase_quality = PhaseDerivativeVariance(im_phase);   

%% Unwrap
colref = startingpoint(2); 
rowref = startingpoint(1);
im_unwrapped(rowref,colref) = im_phase(rowref,colref);                      %Save the unwrapped values
unwrapped_binary(rowref,colref,1) = 1;
%Mark the pixels adjoining the selected point
if im_mask(rowref-1, colref, 1)==1 
    adjoin(rowref-1, colref, 1)=1; 
end       
if im_mask(rowref+1, colref, 1)==1 
    adjoin(rowref+1, colref, 1)=1; 
end
if im_mask(rowref, colref-1, 1)==1 
    adjoin(rowref, colref-1, 1)=1; 
end
if im_mask(rowref, colref+1, 1)==1 
    adjoin(rowref, colref+1, 1)=1; 
end
im_unwrapped = GuidedFloodFill(im_phase, im_unwrapped, ...
    unwrapped_binary, im_phase_quality, adjoin, im_mask);    %Unwrap

end