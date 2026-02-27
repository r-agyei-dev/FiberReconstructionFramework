%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MAIN CODE FOR THE 2D SEGMENTATION OF TOMOGRAMS
% BY RONALD AGYEI
% DATE:
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DESCRIPTION:
% This function performs automated 2D segmentation on chunked tomographic
% image volumes. The workflow includes:
%   1. Parallel pool initialization
%   2. Iterative sharpening + watershed segmentation
%   3. Comparative analysis of segmentation results
%   4. Secondary watershed refinement
%   5. Saving and final merging of chunk results
%
% INPUT ARGUMENTS (varargin):
%   1) tomo_dataset_idx : Dataset identifier string
%   2) num_work         : Number of parallel workers
%   3) mfile_path_file  : Path to main m-file directory
%   4) save_path        : Path containing chunked data
%   5) flip_axes_var    : Axis configuration flag (1=XYZ,2=ZXY,3=YZX)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [varargout] = run2DSegmentationPipeline(varargin)

%% =========================================================================
% INPUT ARGUMENT PARSING
% =========================================================================
tomo_dataset_idx      = varargin{1};   % Dataset index for specimen
num_work              = varargin{2};   % Number of MATLAB workers
mfile_path_file       = varargin{3};   % Path to m-file directory
save_path             = varargin{4};   % Path to chunked dataset
flip_axes_var         = varargin{5};   % Axis configuration selector

%% =========================================================================
% STEP 1: Initiate MATLAB Parallel Pool
% =========================================================================
disp('% STEP 1: Initiate Matlab Parallelization using PARPOOL')

if isempty(gcp('nocreate'))
    myCluster = parcluster('local');
    myCluster.NumWorkers = num_work;
    parpool(myCluster);
end

%% =========================================================================
% STEP 2: Define Output Naming and Locate Chunk Files
% =========================================================================
disp('BEGIN SEGMENTATION XYZ')

% Determine final merged output filename based on axis ordering
if flip_axes_var == 1
    main_variable_name = sprintf('%s%s%s', ...
        'FULL_VOLUME_METADATA_LIBRARY_',tomo_dataset_idx,'_XYZ.mat');
elseif flip_axes_var == 2
    main_variable_name = sprintf('%s%s%s', ...
        'FULL_VOLUME_METADATA_LIBRARY_',tomo_dataset_idx,'_ZXY.mat');
elseif flip_axes_var == 3
    main_variable_name = sprintf('%s%s%s', ...
        'FULL_VOLUME_METADATA_LIBRARY_',tomo_dataset_idx,'_YZX.mat');
end

cd(mfile_path_file);

% Locate chunk files depending on axis configuration
if flip_axes_var == 1
    SPLIT_PATH_var = sprintf('%s%s',save_path,'/XYZ_files');
elseif flip_axes_var == 2
    SPLIT_PATH_var = sprintf('%s%s',save_path,'/ZXY_files');
elseif flip_axes_var == 3
    SPLIT_PATH_var = sprintf('%s%s',save_path,'/YZX_files');
end

patternname = sprintf('%s%s%s%s','Input_',tomo_dataset_idx,'_','*.mat');
filelist    = dir(fullfile(SPLIT_PATH_var,patternname));
num_chunks  = numel(filelist);

%% =========================================================================
% STEP 3: Parameter Space Definition for Segmentation
% =========================================================================
%
% These parameters control:
%   • imsharpen
%   • imadjust
%   • watershed marker generation
%
% All combinations are generated using ndgrid.

% ---- imsharpen parameters ----
amount_thresh = 1.5;     % Sharpen strength (typical range 0–2)
radius        = 2;       % Neighborhood radius for sharpening
threshold     = 0.2;     % Edge threshold

% ---- imadjust parameters ----
im_adj_high   = 0.6:0.1:0.9;  % Upper intensity mapping range
im_adj_low    = 0.2:0.1:0.4;  % Lower intensity mapping range

% ---- watershed parameter ----
Orig_level    = 0.7:0.1:0.8;  % Marker sensitivity

% Generate full parameter grid
A = amount_thresh;
B = radius;
C = threshold;
H = im_adj_high;
L = im_adj_low;
O = Orig_level;

[ca, cb, cc, ch, cl, co] = ndgrid(A,B,C,H,L,O);
d_amt_rad_adj = [ca(:), cb(:), cc(:), ch(:), cl(:), co(:)];

%% =========================================================================
% STEP 4: Process Each Chunk Independently
% =========================================================================
for sp_idx = 1:num_chunks

    % ----- Construct per-chunk output filename -----
    if flip_axes_var == 1
        main_variable_name_pt = sprintf('%s%s%s%d%s', ...
            'FULL_VOLUME_METADATA_LIBRARY_',tomo_dataset_idx,...
            '_XYZ_PART_',sp_idx,'.mat');
    elseif flip_axes_var == 2
        main_variable_name_pt = sprintf('%s%s%s%d%s', ...
            'FULL_VOLUME_METADATA_LIBRARY_',tomo_dataset_idx,...
            '_ZXY_PART_',sp_idx,'.mat');
    elseif flip_axes_var == 3
        main_variable_name_pt = sprintf('%s%s%s%d%s', ...
            'FULL_VOLUME_METADATA_LIBRARY_',tomo_dataset_idx,...
            '_YZX_PART_',sp_idx,'.mat');
    end

    % ----- Load current chunk -----
    image_vol_string = sprintf('%s%s%s%d%s', ...
        'Input_',tomo_dataset_idx,'_',sp_idx,'.mat');

    Input = load(sprintf('%s%s%s',SPLIT_PATH_var,'/',image_vol_string));
    Input = cell2mat(struct2cell(Input));

%% =========================================================================
% STEP 4A: Iterative Sharpening + Watershed Segmentation
% =========================================================================
    aa                 = 1;
    numberOfImageFiles = size(Input,3);
    num_slices         = numberOfImageFiles - 1;
    size_length_2D     = [size(Input,1) size(Input,2)];

    disp('BEGIN SHARPENING & WATERSHED SEGMENTATION PROCEDURE FOR IMAGES')

    [columnsInImage,rowsInImage] = ...
        meshgrid(1:size_length_2D(2),1:size_length_2D(1));

    Input = uint16(Input);
    MATLAB_SHARP_comp = cell(1,numel(aa:aa+num_slices));

    % ----- Loop through slices -----
    for ii = aa:aa+num_slices

        sprintf('%s%d','Sharpen & Segment Image Slice ',ii)

        if ~isempty(nonzeros(Input(:,:,ii)))

            IMAGE_TEMP = Input(:,:,ii);

            % =========================================================
            % Parallel sweep over parameter combinations
            % =========================================================
            parfor kk = 1:size(d_amt_rad_adj,1)

                if kk == 1
                    Adjust = 0;  % No intensity adjustment
                else
                    Adjust = 1;  % Apply intensity adjustment
                end

                [MATLAB_SHARP{kk},~] = ...
                    ITERATIVE_IM_SHARP_WATERSHED_SEG_FUNCTION( ...
                    IMAGE_TEMP,...
                    d_amt_rad_adj(kk,:),...
                    Adjust,...
                    size_length_2D,...
                    columnsInImage,...
                    rowsInImage);
            end

            % Remove empty results
            MATLAB_SHARP( ...
                cell2mat(cellfun(@(x)isempty(x),MATLAB_SHARP,'uni',0)) ...
                ) = [];

            MATLAB_SHARP_comp{ii} = MATLAB_SHARP';
        end
    end

    disp('END SHARPENING & WATERSHED SEGMENTATION PROCEDURE FOR IMAGES')

%% =========================================================================
% STEP 5: Comparative Analysis Across Parameter Results
% =========================================================================
    SLICE_REGIONS_STATS_STRUCTURE = ...
        cell(1,numel(aa:aa+num_slices));

    SLICE_REGIONS_STATS_TABLE = ...
        cell(1,numel(aa:aa+num_slices));

    parfor ii = aa:aa+num_slices

        sprintf('%s%d', ...
            'BEGIN COMPARATIVE ANALYSIS FOR SHARPENED IMAGES',ii)

        if ~isempty(nonzeros(Input(:,:,ii)))

            MATLAB_TEMP = MATLAB_SHARP_comp{ii}{1};

            % Compare all segmentation variants
            for kk = 2:length(MATLAB_SHARP_comp{ii})

                if kk == 2
                    Comp_Bin_Im       = zeros(size(Input(:,:,ii)));
                    Base_AREA_compare = 0;
                    orig_comp         = 1;
                else
                    orig_comp         = 0;
                end

                [MATLAB_TEMP,Comp_Bin_Im,Base_AREA_compare] = ...
                    COMPARATIVE_ANALYSIS_OF_SEG_RESULTS_A( ...
                    MATLAB_TEMP,...
                    MATLAB_SHARP_comp{ii}{kk},...
                    Input(:,:,ii),...
                    Comp_Bin_Im,...
                    orig_comp,...
                    Base_AREA_compare,...
                    columnsInImage,...
                    rowsInImage);
            end

            sprintf('%s%d', ...
                'END COMPARATIVE ANALYSIS FOR SHARPENED IMAGES',ii)

            SLICE_REGIONS_STATS_STRUCTURE{ii} = MATLAB_TEMP;
        end
    end

%% =========================================================================
% STEP 6: Secondary Watershed Refinement
% =========================================================================
    mask_value = 0.8;
    Stats_cell = cell(1,size(Input,3));

    for pp = 1:size(Input,3)

        Image_crop = Input(:,:,pp);

        if ~isempty(nonzeros(Input(:,:,pp)))

            level = graythresh(Image_crop);

            Stats_cell{pp} = ...
                SECONDARY_WATERSHED_FUNCTION( ...
                Image_crop,level,mask_value);

            SLICE_REGIONS_STATS_STRUCTURE{pp} = ...
                COMPARATIVE_ANALYSIS_OF_SEG_RESULTS_B( ...
                SLICE_REGIONS_STATS_STRUCTURE{pp},...
                Stats_cell{pp},...
                size_length_2D);

            SLICE_REGIONS_STATS_TABLE{pp} = ...
                struct2table(SLICE_REGIONS_STATS_STRUCTURE{pp});
        end
    end

%% =========================================================================
% STEP 7: Save Chunk Results
% =========================================================================
    disp('BEGIN SAVING GLOBAL VARIABLES FOR RECONSTRUCTION')

    size_length = [size(Input,1) size(Input,2) ...
                   length(SLICE_REGIONS_STATS_STRUCTURE)];

    save(sprintf('%s',main_variable_name_pt), ...
        'SLICE_REGIONS_STATS_TABLE','-v7.3')

    movefile(main_variable_name_pt,SPLIT_PATH_var)

    disp('END SAVING GLOBAL  VARIABLES FOR RECONSTRUCTION')
end

%% =========================================================================
% FINAL MERGE OF ALL CHUNKS
% =========================================================================
clearvars -except SPLIT_PATH SPLIT_PATH_var tomo_dataset_idx ...
    flip_axes_var num_chunks size_length main_variable_name ...
    save_path mfile_path_file wrong_reconstruction size_length

cd(SPLIT_PATH_var)

if flip_axes_var == 1
    main_variable_name_pre = ...
        sprintf('%s%s%s', ...
        'FULL_VOLUME_METADATA_LIBRARY_',tomo_dataset_idx,'_XYZ_PART_');
elseif flip_axes_var == 2
    main_variable_name_pre = ...
        sprintf('%s%s%s', ...
        'FULL_VOLUME_METADATA_LIBRARY_',tomo_dataset_idx,'_ZXY_PART_');
elseif flip_axes_var == 3
    main_variable_name_pre = ...
        sprintf('%s%s%s', ...
        'FULL_VOLUME_METADATA_LIBRARY_',tomo_dataset_idx,'_YZX_PART_');
end

SLICE_REGIONS_STATS_STRUCTURE_cell = cell(1,num_chunks);

for tt = 1:num_chunks
    load(sprintf('%s%d%s',main_variable_name_pre,tt,'.mat'));
    SLICE_REGIONS_STATS_STRUCTURE_cell{tt} = ...
        SLICE_REGIONS_STATS_TABLE;
    clearvars SLICE_REGIONS_STATS_STRUCTURE_t
end

SLICE_REGIONS_STATS_STRUCTURE = ...
    horzcat(SLICE_REGIONS_STATS_STRUCTURE_cell{:});

size_length = [size_length(1) size_length(2) ...
               size(SLICE_REGIONS_STATS_STRUCTURE,2)];

cd(mfile_path_file);

save(main_variable_name, ...
    'SLICE_REGIONS_STATS_STRUCTURE','size_length','-v7.3');

movefile(main_variable_name,save_path)

varargout{1} = SLICE_REGIONS_STATS_STRUCTURE;

end