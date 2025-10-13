clc; close all; clear;
%% -------------------------- Path Setup --------------------------------------
currentFolder = pwd;
[~, folderName] = fileparts(currentFolder);

% Ensure correct working directory
if ~strcmp(folderName, 'Code5_Phase_retrieval')
    error('Run this script from "Code5_Phase_retrieval". Current folder: %s', folderName);
end

% Add required functions (relative path instead of absolute)
rootFolder = fileparts(currentFolder);
addpath(genpath(fullfile(currentFolder, 'data')));

%% -------------------------- Config ------------------------------------------
% File names for new and previous image_stack
new_syscor_name = 'pr1_syscorr_p2um';   % will be saved/updated
pre_syscor_name = 'PR0_flat';           % baseline syscor
defocus = 0.0;                          % defocus term (um)

% -------------------------- Run Phase Retrieval ------------------------------
phase_retrieval_defocus_Widefield_hex(new_syscor_name, pre_syscor_name, defocus);

disp('Phase retrieval update complete.');
