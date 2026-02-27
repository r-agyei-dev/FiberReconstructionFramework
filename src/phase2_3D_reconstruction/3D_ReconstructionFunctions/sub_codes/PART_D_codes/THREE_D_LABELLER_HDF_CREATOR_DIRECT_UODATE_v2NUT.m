%% THREE_D_LABELLER_HDF_CREATOR_DIRECT_UODATE_v2NUT
% This function generates a 3D labeled volume for fiber datasets, handling
% optional sieving of fibers based on angular thresholds. It outputs the
% labeled volume, fiber heights, linear indices, z-coordinate vectors,
% corrected centroids, out-of-plane angles, and missing fiber indices.
%
% Differences from v2 version:
%   - No split checking or detailed fiber refinement
%   - Simpler version focused on creating labeled volume and HDF export
%
% Inputs:
%   check_mat_center_GLOBAL        - Cell array of linear indices for each fiber
%   size_length                    - Dimensions of the 3D volume [X,Y,Z]
%   lump                           - Flag to combine fibers for height calculation
%   hdfsave                        - Flag to save output as HDF
%   Centroid_MID_cell_LUMP_UPDATE  - Cell array of fiber centroids
%   load_variable                  - Determines voxel ordering: 1: XYZ, 2: ZXY, 3: YZX
%   FINAL_CORRECTED_CENTROID       - Fiber centroid data (already corrected)
%   new_idx                        - Dataset version/index for output
%   save_path                      - Directory for HDF files
%   sieve_num                      - Identifier for HDF indexing
%   grad_out_LUMP_UPDATE            - Out-of-plane angles for fibers
%   sieve_indicator                - Flag to remove fibers above angular threshold
%   angular_thresh                 - Angular threshold for sieving fibers
%
% Outputs:
%   varargout{1}  - 3D labeled volume (TEMP_MAT)
%   varargout{2}  - Fiber height and index matrix (Idx_Height)
%   varargout{3}  - Linear indices of each fiber (FIB_LIN_IDX)
%   varargout{4}  - Z-coordinate vectors of each fiber (z_vec)
%   varargout{5}  - Corrected fiber centroids (FINAL_CORRECTED_CENTROID)
%   varargout{6}  - Out-of-plane angles (grad_out_LUMP_UPDATE)
%   varargout{7}  - Indices of missing fibers (missing_idx)

function [varargout] = THREE_D_LABELLER_HDF_CREATOR_DIRECT_UODATE_v2NUT(varargin)

%% Load input variables
check_mat_center_GLOBAL           = varargin{1};
size_length                       = varargin{2};
lump                              = varargin{3};
hdfsave                           = varargin{4};
Centroid_MID_cell_LUMP_UPDATE     = varargin{5};
load_variable                     = varargin{6};
FINAL_CORRECTED_CENTROID          = varargin{7};
new_idx                           = varargin{8};
save_path                          = varargin{9};
sieve_num                          = varargin{10};
grad_out_LUMP_UPDATE               = varargin{11};
sieve_indicator                     = varargin{12};
angular_thresh                     = varargin{13};

%% Optional sieving of fibers based on angular threshold
if sieve_indicator == 1
    % Remove fibers exceeding the threshold angle
    Centroid_MID_cell_LUMP_UPDATE(grad_out_LUMP_UPDATE > angular_thresh) = [];
    FINAL_CORRECTED_CENTROID(grad_out_LUMP_UPDATE > angular_thresh) = [];
    check_mat_center_GLOBAL(grad_out_LUMP_UPDATE > angular_thresh) = [];
    grad_out_LUMP_UPDATE(grad_out_LUMP_UPDATE > angular_thresh) = [];
    
    nom_bin_prefix  = 'REFINED_FIBERS_SIEVED';
    nom_file_prefix = 'REFINED_FIBERS_SIEVED';
else
    nom_bin_prefix  = 'REFINED_FIBERS_UNSIEVED';
    nom_file_prefix = 'REFINED_FIBERS_UNSIEVED';
end

%% Remove empty fibers remaining after sieving
empty_mask = cell2mat(cellfun(@(x) isempty(x), check_mat_center_GLOBAL, 'uni', 0));
Centroid_MID_cell_LUMP_UPDATE(empty_mask) = [];
FINAL_CORRECTED_CENTROID(empty_mask) = [];
grad_out_LUMP_UPDATE(empty_mask) = [];
check_mat_center_GLOBAL(empty_mask) = [];

%% Compute fiber height index mapping
if lump == 1
    % Combine all fibers into a single vector for labeling
    check_mat_center_GLOBAL_VEC = [check_mat_center_GLOBAL{:}];
    Idx_Height = [];
else
    % Keep fibers as separate cells
    check_mat_center_GLOBAL_VEC = check_mat_center_GLOBAL;
    Idx_Height = cell2mat(cellfun(@(x) size(x,2), Centroid_MID_cell_LUMP_UPDATE, 'uni', 0))';
end
Idx_Height(:,2) = 1:size(Idx_Height(:,1),1);
Idx_Height = fliplr(Idx_Height);

%% Initialize 3D labeled volume and fiber arrays
TEMP_MAT = zeros([max(size_length), max(size_length), min(size_length)]);
FIB_LIN_IDX = cell(1, length(check_mat_center_GLOBAL_VEC));
z_vec = cell(1, length(check_mat_center_GLOBAL_VEC));

%% Assign fiber indices to 3D volume
for ii = 1:length(check_mat_center_GLOBAL_VEC)
    [x, y, z] = ind2sub(size_length, check_mat_center_GLOBAL_VEC{ii});
    z_vec{ii} = z(1):z(end);
    xyz = [x, y, z];
    
    % Permute axes if volume ordering is non-standard
    if load_variable == 2  % ZXY
        perm_idx = [2 3 1];
    elseif load_variable == 3  % YZX
        perm_idx = [3 1 2];
    else
        perm_idx = [1 2 3];
    end
    
    if load_variable > 1
        size_length_p = size_length(perm_idx);
        xyz = xyz(:, perm_idx);
    else
        size_length_p = size_length;
    end
    
    % Convert to linear indices and fill 3D volume
    FIB_LIN_IDX{ii} = sub2ind(size_length_p, xyz(:,1), xyz(:,2), xyz(:,3));
    TEMP_MAT(FIB_LIN_IDX{ii}) = ii;
end

%% Identify missing fibers (indices not represented in volume)
missing_idx = find(~ismember(1:size(check_mat_center_GLOBAL_VEC,2), nonzeros(unique(TEMP_MAT(:)))));

%% Generate HDF file for the labeled volume
[~] = MAIN_HDF_CREATOR(TEMP_MAT, nom_bin_prefix, nom_file_prefix, sieve_num, new_idx, hdfsave, load_variable, save_path);

%% Assign outputs
varargout{1} = TEMP_MAT;
varargout{2} = Idx_Height;
varargout{3} = FIB_LIN_IDX;
varargout{4} = z_vec;
varargout{5} = FINAL_CORRECTED_CENTROID;
varargout{6} = grad_out_LUMP_UPDATE;
varargout{7} = missing_idx;

end