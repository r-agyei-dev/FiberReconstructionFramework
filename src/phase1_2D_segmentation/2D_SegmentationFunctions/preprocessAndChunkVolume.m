function preprocessAndChunkVolume( ...
    images_file_directory, ...
    img_vol_chunks_directory, ...
    mat_file_directory, ...
    grey_volume_name, ...
    tomo_dataset_idx, ...
    flip_axes_var, ...
    number_per_chunck, ...
    Data_Crop, ...
    crop_var, ...
    tiff_or_mat)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ALGORITHM FOR PRE-PROCESSING (Function Version)
% FINAL VERSION UPDATE: 10/04/2018
% BY: Ronald F Agyei
%
% Converts large grayscale volumes into chunked .mat files
% to mitigate MATLAB memory limitations.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

mfile_path_file = pwd;

%% ------------------------------------------------------------------------
% Ensure output directories exist
% -------------------------------------------------------------------------
if ~exist(mat_file_directory,'dir')
    mkdir(mat_file_directory)
end

if ~exist(img_vol_chunks_directory,'dir')
    mkdir(img_vol_chunks_directory)
end

%% ------------------------------------------------------------------------
% Load image volume
% -------------------------------------------------------------------------
fileFolder = images_file_directory;

if tiff_or_mat == 1
    image_type_ext = '*.tif';
    Input = Tiff_Fast_Loader(fileFolder, image_type_ext);
else
    s  = load(fullfile(mat_file_directory,grey_volume_name));
    fn = fieldnames(s);
    Input = s.(fn{1});
end

%% ------------------------------------------------------------------------
% Optional cropping
% -------------------------------------------------------------------------
if crop_var == 1
    Input = Input( ...
        Data_Crop.ycenter-Data_Crop.radius : Data_Crop.ycenter+Data_Crop.radius-1, ...
        Data_Crop.xcenter-Data_Crop.radius : Data_Crop.xcenter+Data_Crop.radius-1, ...
        :);
end

%% ------------------------------------------------------------------------
% Save full grey volume
% -------------------------------------------------------------------------
grey_name = sprintf('Grey_Vol_%s.mat', tomo_dataset_idx);
save(grey_name, 'Input', '-v7.3');
movefile(grey_name, mat_file_directory);

%% ------------------------------------------------------------------------
% Loop over flip axes
% -------------------------------------------------------------------------
for xx = 1:numel(flip_axes_var)

    % -------- Determine orientation --------
    switch flip_axes_var(xx)

        case 1
            SPLIT_PATH = fullfile(img_vol_chunks_directory,'XYZ_files');

        case 2
            Input = permute(Input,[3 1 2]); % xyz -> zxy
            SPLIT_PATH = fullfile(img_vol_chunks_directory,'ZXY_files');

        case 3
            Input = permute(Input,[2 3 1]); % xyz -> yzx
            SPLIT_PATH = fullfile(img_vol_chunks_directory,'YZX_files');
    end

    if ~exist(SPLIT_PATH,'dir')
        mkdir(SPLIT_PATH)
    end

    cd(mfile_path_file)

    %% --------------------------------------------------------------------
    % Chunk the volume
    % ---------------------------------------------------------------------
    slice_depth = size(Input,3);
    slice_idx   = 1:number_per_chunck:slice_depth;

    for uu = 1:numel(slice_idx)

        if uu == 1
            slice_tmp = 1:slice_idx(uu+1)-1;
        elseif uu < numel(slice_idx)
            slice_tmp = slice_idx(uu):slice_idx(uu+1)-1;
        else
            slice_tmp = slice_idx(uu):slice_depth;
        end

        chunk_data = Input(:,:,slice_tmp);

        var_name = sprintf('Input_%s_%d', tomo_dataset_idx, uu);
        file_name = sprintf('%s.mat', var_name);

        % Save without eval
        S.(var_name) = chunk_data;
        save(file_name,'-struct','S');
        S = struct(); % reset

        movefile(file_name, SPLIT_PATH);
    end

    %% --------------------------------------------------------------------
    % Restore original orientation
    % ---------------------------------------------------------------------
    switch flip_axes_var(xx)
        case 2
            Input = permute(Input,[2 3 1]); % back to xyz
        case 3
            Input = permute(Input,[3 1 2]); % back to xyz
    end

end

end