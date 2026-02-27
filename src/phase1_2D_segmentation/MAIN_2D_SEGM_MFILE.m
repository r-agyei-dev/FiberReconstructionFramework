%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 2D SEGMENTATION PIPELINE FOR SERIAL TOMOGRAPHIC SECTIONS
%
% Final Version: 02/22/2026
% Author: Ronald F. Agyei
%
% DESCRIPTION:
% This script performs automated 2D segmentation on serial slices extracted
% from a 3D tomographic image volume. Segmentation is executed independently
% along specified principal axes (x, y, z) to support robust downstream
% 3D reconstruction.
%
% PIPELINE OVERVIEW:
%   1. Define dataset and computational parameters.
%   2. Preprocess and chunk the grayscale volume.
%   3. Initialize parallel segmentation.
%   4. Perform axis-wise 2D segmentation.
%   5. Save segmented outputs to the results directory.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars; close all; clc;

%% =========================================================================
% USER PARAMETERS
% =========================================================================

tomo_dataset_idx          = '0002';          % Dataset identifier
num_workers               = 12;              % Parallel workers
flip_axes_var             = 1:3;             % 1=x, 2=y, 3=z
mfile_path_file           = pwd;

% -------------------- Paths --------------------
images_file_directory     = '/Users/ronaldagyei/Desktop/PhD_Codes/results';
img_vol_chunks_directory  = '/Users/ronaldagyei/Desktop/PhD_Codes/results/SPLIT_FILES';
mat_file_directory        = '/Users/ronaldagyei/Desktop/PhD_Codes/mat_files';
save_path                 = img_vol_chunks_directory;

% -------------------- Volume Settings --------------------
grey_volume_name          = 'New_Grey_Data_S2.mat';
number_per_chunck         = 30;
tiff_or_mat               = 2;               % 1=tiff, 2=mat
crop_var                  = 0;               % 1=enable cropping

% -------------------- Crop Parameters --------------------
Data_Crop.xMeshgrid       = 2560;
Data_Crop.yMeshgrid       = 2560;
Data_Crop.radius          = 946;
Data_Crop.xcenter         = 1328;
Data_Crop.ycenter         = 1283;

%% =========================================================================
% PREPROCESS AND CHUNK VOLUME
% =========================================================================

preprocessAndChunkVolume( ...
    images_file_directory, ...
    img_vol_chunks_directory, ...
    mat_file_directory, ...
    grey_volume_name, ...
    tomo_dataset_idx, ...
    flip_axes_var, ...
    number_per_chunck, ...
    Data_Crop, ...
    crop_var, ...
    tiff_or_mat);

%% =========================================================================
% ENSURE OUTPUT DIRECTORY EXISTS
% =========================================================================

if ~exist(save_path,'dir')
    mkdir(save_path);
end

%% =========================================================================
% RUN 2D SEGMENTATION ALONG SPECIFIED AXES
% =========================================================================

fprintf('\nStarting 2D segmentation...\n');

% tic
% SEGMENTATION_RESULTS = run2DSegmentationPipeline(...
%                                                  tomo_dataset_idx,...
%                                                  save_path,...
%                                                  img_vol_chunks_directory,...
%                                                  num_workers);
% 
% toc

fprintf('Segmentation pipeline completed.\n');

tic

for ii = 1:numel(flip_axes_var)

    current_axis = flip_axes_var(ii);
    fprintf('Processing axis %d...\n', current_axis);

    run2DSegmentationPipeline( ...
        tomo_dataset_idx, ...
        num_workers, ...
        mfile_path_file,...
        save_path, ...
        current_axis);

end

toc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%