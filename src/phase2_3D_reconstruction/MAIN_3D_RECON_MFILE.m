%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% ALGORITHM FOR 3D RECONSTRUCTION CODE
% FINAL VERSION UPDATE: 9/20/2018
% BY: Ronald F. Agyei
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% This script performs 3D fiber reconstruction from 2D segmented elliptical 
% features extracted from greyscale tomograms. It handles orientation, 
% stacking, correction of reconstruction errors, and saving the refined 3D microstructure.

clearvars; clc; close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% DESCRIPTION OF VARIABLES
% Mfile_Directory          : Directory where the MATLAB code resides
% tomo_dataset_idx         : Index (string) of tomogram dataset for the specimen
% num_work                 : Number of MATLAB parallel workers (not explicitly used here)
% flip_axes_var            : Vector specifying axes along which segmentation occurs (1=x,2=y,3=z)
% mfile_path_file          : Directory path for the current MATLAB file
% specimen_name            : Prefix/name of the specimen
% save_path                : Directory for storing reconstructed images
% angular_thresh           : Threshold angle for filtering spurious fiber inclinations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loop index for reconstructed images (currently single loop)
loop_vector = 1;

% Set directories and parameters
Mfile_Directory   = pwd;                      % Current folder as main directory
flip_axes_var     = 1:3;                      % Reconstruct along x, y, z axes
angular_thresh    = 60;                       % Max angle beyond which fibers are rejected
tomo_dataset_idx  = '0002';                   % Tomogram dataset ID
specimen_name     = 'specimen3_gfrp_';       % Specimen name prefix

% Directory for saving .mat files and sub-code paths
split_file_directory = '/Users/ronaldagyei/Desktop/PhD_Codes/results/SPLIT_FILES/';
mat_file_directory = '/Users/ronaldagyei/Desktop/PhD_Codes/mat_files';
pathAcodes         = '/sub_codes/PART_A_codes';  % Crude microstructure generation
pathBcodes         = '/sub_codes/PART_B_codes';  % 3D fiber reconstruction & stitching
pathDcodes         = '/sub_codes/PART_D_codes';  % Saving refined microstructure

% Cropping and slicing parameters
crop_vol     = 1;    % Flag to crop volume
start_slice  = 1;    % Starting slice for reconstruction

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% LOOP OVER DATASETS
for pp = 1
    new_idx = loop_vector(pp);  % Current dataset index
    
    % Loop over three axes: x=1, y=2, z=3
    for kk = 1:3
        load_variable = kk;   % Determines current axis
        
        % Directory to save reconstructed files for current specimen
        % fileFolder_1 = sprintf('%s%s%s%s', split_file_directory, '/RECON_FILES/', specimen_name, tomo_dataset_idx);
        % mkdir(fileFolder_1)
        
        % Directories and file prefixes for metadata
        segmat_prenom = '/FULL_VOLUME_METADATA_LIBRARY_';
        
        % Load the correct orientation of 2D segmentation metadata based on axis
        if load_variable == 1
            load(sprintf('%s%s%s%s', split_file_directory, segmat_prenom, tomo_dataset_idx, '_XYZ.mat'))
            save_path = sprintf('%s%s', split_file_directory, 'XYZ_files');
        elseif load_variable == 2 
            load(sprintf('%s%s%s%s', split_file_directory, segmat_prenom, tomo_dataset_idx, '_ZXY.mat'))
            save_path = sprintf('%s%s', split_file_directory, 'ZXY_files');
        elseif load_variable == 3
            load(sprintf('%s%s%s%s', split_file_directory, segmat_prenom, tomo_dataset_idx, '_YZX.mat'))
            save_path = sprintf('%s%s', split_file_directory, 'YZX_files');
        end
        
        % % Create extra folder for intermediate files
        % save_path_extra_file = sprintf('%s%s%d', save_path, '/Extra_Files');
        % mkdir(save_path)
        % mkdir(save_path_extra_file)
        
        % Return to main directory
        cd(Mfile_Directory)
        
        %% STEP 1: Binarize RAW microstructure
        disp('Create the Binarized version of RAW microstructure')
        
        % Load greyscale tomogram data
        TT = load(sprintf('%s%s%s', mat_file_directory, '/New_Grey_Data_S2.mat'));
        names = fieldnames(TT);
        Image_vol = TT.(names{1});    % Extract volume data
        clearvars TT
        
        orig_size_length = size(Image_vol);  % Original size for later use
        
        % Permute volume according to current axis
        if load_variable == 2
            Image_vol = permute(Image_vol, [3 1 2]);  % ZXY orientation
        elseif load_variable == 3
            Image_vol = permute(Image_vol, [2 3 1]);  % YZX orientation
        end
        
        size_length = size(Image_vol);  % Size after permutation
        
        % Generate binarized volume from greyscale
        nom_bin_prefix = sprintf('%s%d', 'Actual_Binary_Image_vol_', load_variable);
        [~, Pixels_NON_BIN] = BINARIZED_GREYSCALE_VOL_FUNCTION(Image_vol, nom_bin_prefix, save_path, tomo_dataset_idx, load_variable);
        
        disp('Clean the STRUCTURE ARRAY before elliptical stacking')
        
        % Convert slice statistics tables to structure array
        for tt = 1:size(SLICE_REGIONS_STATS_STRUCTURE, 2)-1
            SLICE_REGIONS_STATS_STRUCTURE{tt} = table2struct(SLICE_REGIONS_STATS_STRUCTURE{tt}); 
        end
        
        % Extract centroids for each slice
        Centroid_cell = cell(1, length(SLICE_REGIONS_STATS_STRUCTURE));
        for ii = 1:length(SLICE_REGIONS_STATS_STRUCTURE)-1
            Centroid_cell{ii} = cat(1, (SLICE_REGIONS_STATS_STRUCTURE{ii}(:).Centroid));   
        end
        clearvars Image_vol
        
        %% STEP 2: Generate crude microstructure by stacking ellipses
        if crop_vol == 1 && load_variable == 1
            numberOfImageFiles = size(SLICE_REGIONS_STATS_STRUCTURE,2)-1; 
        else
            numberOfImageFiles = size_length(3)-1;
        end
        iter_limit = numberOfImageFiles - 1;  % Limit for stacking
        
        % Add PART A codes path
        addpath(sprintf('%s%s', Mfile_Directory, pathAcodes))
        
        % Generate crude 3D fiber microstructure
        [done, properties, idx_crude_fibers, cat_size_const_ellipses] = ...
            CRUDE_MICROSTRUCTURE_GENERATION_FUNCTION(SLICE_REGIONS_STATS_STRUCTURE, Centroid_cell, ...
                                                      numberOfImageFiles, size_length, iter_limit);
        
        %% STEP 3: Correct erroneous 3D reconstructions (e.g., conjoined fibers)
        addpath(sprintf('%s%s', Mfile_Directory, pathBcodes))
        ORIGINAL_CRUDE_BIN_MAT = 0; % Placeholder/check
        
        [Linear_Index_Final, Pixel_Idx_Per_Slice_MAIN, ...
         FINAL_CORRECTED_CENTROID, IDX_to_keep, OUT_PLANE_ANGLE, ...
         height_variable_MAIN, fib_diam_stitched_MAIN_VARIABLE] = ...
            FIBER_RECONSTRUCTION_CORRECTION_FUNCTION(properties, iter_limit, size_length, load_variable, ...
                                                     angular_thresh, save_path, new_idx);
        
        disp('Saving before Section D')
        
        %% STEP 3B: Perform fiber splitting to separate connected fibers
        clearvars cat_size_const_ellipses Centroid_cell idx_crude_fibers IDX_to_keep
        clearvars OUT_PLANE_ANGLE properties SLICE_REGIONS_STATS_STRUCTURE
        
        [split_array_slice, split_array_consec_int] = Fiber_Stitch_analysis_CP_SIII(Linear_Index_Final, height_variable_MAIN, ...
                                                                                     Pixels_NON_BIN, size_length);
                                                                                 
        % Create temporary copies for splitting
        Linear_Index_Final_TMP       = Linear_Index_Final;
        FINAL_CORRECTED_CENTROID_TMP = FINAL_CORRECTED_CENTROID;
        Pixel_Idx_Per_Slice_MAIN_TMP = Pixel_Idx_Per_Slice_MAIN;
        height_variable_MAIN_TMP     = height_variable_MAIN;
        
        disp('Perform the split here')
        
        % Find fibers that require splitting
        split_fibers = find(cell2mat(cellfun(@(x)~isempty(x), split_array_consec_int, 'uni', 0)));
        
        if ~isempty(split_fibers)
            for i = 1:length(split_fibers)
                idx_tmp = split_fibers(i);
                
                % Split Linear Index, Centroid, Pixel Index, and Height variables
                Linear_Index_Final_TMP{idx_tmp} = cellfun(@(x)Linear_Index_Final{idx_tmp}(ismember(height_variable_MAIN{idx_tmp}, x)), split_array_slice{idx_tmp}, 'uni', 0);
                FINAL_CORRECTED_CENTROID_TMP{idx_tmp} = cellfun(@(x)FINAL_CORRECTED_CENTROID{idx_tmp}(ismember(FINAL_CORRECTED_CENTROID{idx_tmp}(:,3), x), :), split_array_slice{idx_tmp}, 'uni', 0);
                Pixel_Idx_Per_Slice_MAIN_TMP{idx_tmp} = cellfun(@(x)Pixel_Idx_Per_Slice_MAIN{idx_tmp}(ismember(1:size(Pixel_Idx_Per_Slice_MAIN{idx_tmp}, 2), x)), split_array_consec_int{idx_tmp}, 'uni', 0);
                height_variable_MAIN_TMP{idx_tmp} = cellfun(@(x)height_variable_MAIN{idx_tmp}(ismember(height_variable_MAIN{idx_tmp}, x)), split_array_slice{idx_tmp}, 'uni', 0);
            end
        end
        
        clearvars Linear_Index_Final FINAL_CORRECTED_CENTROID Pixel_Idx_Per_Slice_MAIN height_variable_MAIN
        
        % Concatenate split arrays
        split_array_slice_TMP = split_array_slice;
        split_array_slice_TMP(cell2mat(cellfun(@(x)isempty(x), split_array_slice, 'uni', 0))) = {1};
        Pixel_Idx_Per_Slice_MAIN_TMP_CAT = cell(1,3);
        
        j = 1;
        for i = 1:numel(Linear_Index_Final_TMP)
            if ~iscell(Linear_Index_Final_TMP{i})
                Pixel_Idx_Per_Slice_MAIN_TMP_CAT{j} = Pixel_Idx_Per_Slice_MAIN_TMP{i};
                j = j + 1;
            elseif iscell(Linear_Index_Final_TMP{i})
                for mm = 1:numel(Linear_Index_Final_TMP{i})
                    Pixel_Idx_Per_Slice_MAIN_TMP_CAT{j} = Linear_Index_Final_TMP{i}{mm};
                    j = j + 1;
                end
            end
        end
        clearvars Pixel_Idx_Per_Slice_MAIN_TMP
        
        % Concatenate final structures
        if any(cell2mat(cellfun(@(x)iscell(x), Linear_Index_Final_TMP, 'uni', 0)))
            Linear_Index_Final_TMP_CAT        = [Linear_Index_Final_TMP{:}];
            FINAL_CORRECTED_CENTROID_TMP_CAT = [FINAL_CORRECTED_CENTROID_TMP{:}];
            height_variable_MAIN_TMP_CAT     = [height_variable_MAIN_TMP{:}];
        else
            Linear_Index_Final_TMP_CAT        = Linear_Index_Final_TMP;
            FINAL_CORRECTED_CENTROID_TMP_CAT = FINAL_CORRECTED_CENTROID_TMP;
            height_variable_MAIN_TMP_CAT     = height_variable_MAIN_TMP;
        end
        
        clearvars Linear_Index_Final_TMP FINAL_CORRECTED_CENTROID_TMP height_variable_MAIN_TMP
        
        %% STEP 4: Generate HDF files for refined microstructure
        sieve_num = 3;                  % Number of sieve levels (not elaborated here)
        sieve_indicator_vec = 1;        % 1=sieved, 2=unsieved
        addpath(sprintf('%s%s', Mfile_Directory, pathDcodes))
        
        for ii = 1:numel(sieve_indicator_vec)
            sieve_indicator = sieve_indicator_vec(ii);
            
            % Generate refined 3D microstructure and save
            [stat_report_D, missing_idx, TEMP_MAT] = REFINED_MICROSTRUCTURE_FUNCTION(Linear_Index_Final_TMP_CAT, orig_size_length, size_length, ...
                                                                                    Pixel_Idx_Per_Slice_MAIN_TMP_CAT, load_variable, ...
                                                                                    FINAL_CORRECTED_CENTROID_TMP_CAT, tomo_dataset_idx, save_path, ...
                                                                                    height_variable_MAIN_TMP_CAT, start_slice);
        end
    end
end