% =============================================================================
% Description:
% =============================================================================
% This function binarizes a 3D greyscale image volume and optionally saves it 
% as an HDF5/XDMF volume for visualization.
%
% INPUT VARIABLES:
% ----------------
% Image_vol        : 3D greyscale image volume
% nom_bin_prefix   : Prefix for output H5/XDMF files
% save_path        : Directory to save output files
% tomo_dataset_idx : Dataset index (used in file naming)
% load_variable    : Directional reconstruction indicator
%                    1 = XYZ, 2 = ZXY, 3 = YZX
%
% OUTPUT VARIABLES:
% -----------------
% TEMP_MAT        : 3D binary image volume
% Pixels_NON_BIN  : Logical array indicating background voxels (0 in binary)
%
% NOTE:
% - To save time, HDF5/XDMF output is only created when load_variable == 1
% - The function automatically handles axis permutation for directional reconstruction
% =============================================================================

function [varargout] = BINARIZED_GREYSCALE_VOL_FUNCTION(varargin)

% --------------------------
% Extract input variables
% --------------------------
Image_vol        = varargin{1};
nom_bin_prefix   = varargin{2};
save_path        = varargin{3};
tomo_dataset_idx = varargin{4};
load_variable    = varargin{5};

% Preallocate binary volume
TEMP_MAT = zeros(size(Image_vol));

% Convert to 16-bit for processing
Image_vol = uint16(Image_vol);

disp('Binarizing Image Volume')

% --------------------------
% Step 1: Binarize slice by slice
% --------------------------
for ii = 1:size(Image_vol,3)  
    % Extract current slice
    Image_vol_tmp = Image_vol(:,:,ii);
    
    % Remove zero (background) pixels for threshold calculation
    Image_vol_tmp = Image_vol_tmp(Image_vol_tmp>0);
    
    % Apply Otsu threshold to binarize
    TEMP_MAT(:,:,ii) = imbinarize(Image_vol(:,:,ii), graythresh(Image_vol_tmp)); 
end   

% --------------------------
% Step 2: Identify background voxels
% --------------------------
disp('Finding the Linear Index of NON_BIN areas in the large volume')

% Logical array where 1 = background, 0 = fiber
Pixels_NON_BIN = TEMP_MAT == 0;

% Convert to 8-bit for storage
TEMP_MAT = uint8(TEMP_MAT);

% --------------------------
% Step 3: Permute axes for directional reconstruction
% --------------------------
if load_variable == 3
    TEMP_MAT = permute(TEMP_MAT,[3 1 2]);   % Convert to YZX orientation
elseif load_variable == 2
    TEMP_MAT = permute(TEMP_MAT,[2 3 1]);   % Convert to ZXY orientation
end

% --------------------------
% Step 4: Save HDF5 and XDMF volume (only for load_variable == 1)
% --------------------------
if load_variable == 1
    
    disp('Creating the H5 volume')
    
    % File naming
    nom_bin  = sprintf('%s%s%s', nom_bin_prefix, tomo_dataset_idx, '_.h5');
    nom_file = sprintf('%s%s%s', nom_bin_prefix, tomo_dataset_idx, '_.xdmf');
    
    % Remove previous file if it exists
    delete(nom_bin)
    
    % Determine chunk size for HDF5 (25% of volume size for faster access)
    chunckSize = round(0.25 * size(Image_vol)); 

    % Create HDF5 dataset
    h5create(nom_bin, '/CELL_DATA/FIBER_RECON', size(TEMP_MAT), ...
             'Datatype', 'uint8', 'ChunkSize', chunckSize, 'Deflate', 6);                

    % Dimensions
    [nx, ny, nz] = size(TEMP_MAT);
    dz = chunckSize(3);

    % Determine z-slices to write in chunks
    zStart = 1 : dz : nz;
    zEnd   = dz : dz : nz;

    if zStart(end) == nz
        zStart = zStart(1:end-1);
    end
    if zEnd(end) ~= nz
        zEnd = [zEnd, nz];
    end

    % --------------------------
    % Write volume in chunks
    % --------------------------
    for i = 1:numel(zStart)
        dz = zEnd(i) - zStart(i) + 1;
        if dz > 0
            start = [1 1 zStart(i)];              % Starting indices for HDF5 write
            count = [nx ny dz];                   % Size of chunk to write
            data  = TEMP_MAT(:,:,zStart(i):zEnd(i)); % Chunk data
            h5write(nom_bin, '/CELL_DATA/FIBER_RECON', data, start, count);
        else
            disp('dz == 0, cannot write data')
            break
        end
    end

    % --------------------------
    % Generate XDMF for visualization
    % --------------------------
    h5att_name{1} = sprintf('%s%s%s', nom_bin_prefix, '_xyz_', tomo_dataset_idx);
    h5string{1}   = sprintf('%s%s', nom_bin, ':/CELL_DATA/FIBER_RECON');

    file_name = nom_file;
    z_depth   = size(TEMP_MAT,3);
    x_breadth = size(TEMP_MAT,2);
    y_breadth = size(TEMP_MAT,1);

    disp('Creating the XDMF file from the H5 file')                           
    XDMF_TEXT_GENERATOR(file_name, h5string, h5att_name, z_depth, y_breadth, x_breadth);

    % Move H5 and XDMF files to save directory
    disp('Moving the H5 and XDMF files')                
    movefile(nom_bin, save_path)
    movefile(nom_file, save_path)  
end

% --------------------------
% Step 5: Return outputs
% --------------------------
varargout{1} = TEMP_MAT;           % 3D binary volume
varargout{end+1} = Pixels_NON_BIN; % Logical array of background voxels