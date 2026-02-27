%% THREE_D_LABELLER_HDF_CREATOR_DIRECT_UODATE_v2
% This function processes 3D segmented fiber data to produce a labeled
% volume, compute fiber angles, handle missing fibers, and generate HDF
% files for further analysis. It optionally sieves fibers based on
% out-of-plane angles and refines fragmented fibers before final labeling.
%
% Inputs:
%   check_mat_center_GLOBAL        - Cell array of linear indices for all fibers
%   size_length                    - Volume dimensions [X,Y,Z]
%   lump                           - Flag for whether fibers should be lumped together
%   hdfsave                        - Flag to save HDF files
%   Centroid_MID_cell_LUMP_UPDATE  - Cell array of fiber centroid points per slice
%   load_variable                  - Determines voxel ordering: 1: XYZ, 2: ZXY, 3: YZX
%   FINAL_CORRECTED_CENTROID       - Cell array of corrected fiber centroids
%   new_idx                        - Identifier for dataset versioning
%   save_path                      - Path to save generated HDF files
%   sieve_num                      - Used for indexing HDF outputs
%   grad_out_LUMP_UPDATE            - Out-of-plane angles for fibers
%   sieve_indicator                - Flag to sieve fibers based on angular threshold
%   angular_thresh                 - Threshold angle for sieving
%
% Outputs:
%   varargout{1}  - 3D labeled volume (TEMP_MAT)
%   varargout{2}  - Index and height mapping for fibers (Idx_Height)
%   varargout{3}  - Linear indices of fibers after processing (FIB_LIN_IDX)
%   varargout{4}  - Z-coordinate vectors for each fiber (z_vec)
%   varargout{5}  - Corrected fiber centroids (FINAL_CORRECTED_CENTROID)
%   varargout{6}  - Updated out-of-plane angles (grad_out_LUMP_UPDATE)
%   varargout{7}  - Indices of missing fibers (missing_idx)
%   varargout{8}  - XYZ coordinates of fiber voxels (FIB_PIX_CORD)

function [varargout] = THREE_D_LABELLER_HDF_CREATOR_DIRECT_UODATE_v2(varargin)

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

%% Sieve fibers based on angular threshold if required
if sieve_indicator == 1
    % Remove fibers exceeding angular threshold
    Centroid_MID_cell_LUMP_UPDATE(grad_out_LUMP_UPDATE>angular_thresh) = [];
    FINAL_CORRECTED_CENTROID(grad_out_LUMP_UPDATE>angular_thresh) = [];
    check_mat_center_GLOBAL(grad_out_LUMP_UPDATE>angular_thresh) = [];
    grad_out_LUMP_UPDATE(grad_out_LUMP_UPDATE>angular_thresh) = [];
    
    nom_bin_prefix  = 'REFINED_FIBERS_SIEVED';
    nom_file_prefix = 'REFINED_FIBERS_SIEVED';
else
    nom_bin_prefix  = 'REFINED_FIBERS_UNSIEVED';
    nom_file_prefix = 'REFINED_FIBERS_UNSIEVED';
end

%% Remove any remaining empty fibers
empty_mask = cell2mat(cellfun(@(x) isempty(x), check_mat_center_GLOBAL, 'uni', 0));
Centroid_MID_cell_LUMP_UPDATE(empty_mask) = [];
FINAL_CORRECTED_CENTROID(empty_mask) = [];
grad_out_LUMP_UPDATE(empty_mask) = [];
check_mat_center_GLOBAL(empty_mask) = [];

%% Refine fibers for final labeling
if lump == 1
    check_mat_center_GLOBAL_VEC = [check_mat_center_GLOBAL{:}];  % Combine into single vector
else
    check_mat_center_GLOBAL_VEC = check_mat_center_GLOBAL;       % Keep as cell array
end

%% Detect and correct splits in fibers
[check_mat_center_GLOBAL_VEC, FINAL_CORRECTED_CENTROID] = ...
    REFINED_FIBER_SPLIT_CHECKER(Centroid_MID_cell_LUMP_UPDATE, check_mat_center_GLOBAL_VEC, FINAL_CORRECTED_CENTROID, size_length);

%% Remove empty cells after split correction
missing_idx_1 = cell2mat(cellfun(@(x) isempty(x), check_mat_center_GLOBAL_VEC, 'uni', 0));
missing_idx_2 = cell2mat(cellfun(@(x) isempty(x), FINAL_CORRECTED_CENTROID, 'uni', 0));
check_mat_center_GLOBAL_VEC(missing_idx_1 | missing_idx_2) = [];
FINAL_CORRECTED_CENTROID(missing_idx_1 | missing_idx_2) = [];

%% Compute fiber angles (in-plane and out-of-plane)
single_temp = cell2mat(cellfun(@(x) size(x,1), FINAL_CORRECTED_CENTROID, 'uni', 0)) == 1;
in_plane_cell = cellfun(@(x) atand((x(1,2)-x(end,2)) / (x(end,1)-x(1,1))), FINAL_CORRECTED_CENTROID, 'uni', 0);

% Set in-plane angle to 0 for single-point fibers
[in_plane_cell{single_temp}] = deal(0);

% Handle fibers with overlapping start/end points
overlap_mask = cell2mat(cellfun(@(x) norm(x(1,1:2)-x(end,1:2)), FINAL_CORRECTED_CENTROID, 'uni', 0)) == 0;
[in_plane_cell{overlap_mask}] = deal(0);

% Compute out-of-plane angles
out_plane_angle = cellfun(@(x,p) 90 - acosd((x(end,1)-x(1,1)) / (norm(x(end,:)-x(1,:)) * cosd(p))), ...
                          FINAL_CORRECTED_CENTROID, in_plane_cell, 'uni', 0);

% Replace NaNs with 0
[out_plane_angle{cell2mat(cellfun(@(x) isnan(x), out_plane_angle, 'uni', 0))}] = deal(0);
out_plane_angle = cell2mat(out_plane_angle);

%% Update out-of-plane angles and fiber height indices
grad_out_LUMP_UPDATE = out_plane_angle;
Idx_Height = cell2mat(cellfun(@(x) size(x,1), FINAL_CORRECTED_CENTROID, 'uni', 0))';
Idx_Height(:,2) = 1:size(Idx_Height(:,1),1);
Idx_Height = fliplr(Idx_Height);

%% Create 3D labeled volume and store fiber coordinates
TEMP_MAT = zeros([max(size_length), max(size_length), min(size_length)]);
FIB_LIN_IDX  = cell(1, length(check_mat_center_GLOBAL_VEC));
z_vec        = cell(1, length(check_mat_center_GLOBAL_VEC));
FIB_PIX_CORD = cell(1, length(check_mat_center_GLOBAL_VEC));

for ii = 1:length(check_mat_center_GLOBAL_VEC)
    [x, y, z] = ind2sub(size_length, check_mat_center_GLOBAL_VEC{ii});
    z_vec{ii} = z(1):z(end);
    xyz = [x, y, z];
    
    % Permute axes if required based on load_variable
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
    
    % Store voxel coordinates for each fiber
    FIB_PIX_CORD{ii} = xyz;
end

%% Identify missing fibers
missing_idx = find(~ismember(1:size(check_mat_center_GLOBAL_VEC,2), nonzeros(unique(TEMP_MAT(:)))));

%% Generate HDF files
[~] = MAIN_HDF_CREATOR(TEMP_MAT, nom_bin_prefix, nom_file_prefix, sieve_num, new_idx, hdfsave, load_variable, save_path);

%% Assign outputs
varargout{1} = TEMP_MAT;
varargout{2} = Idx_Height;
varargout{3} = FIB_LIN_IDX;
varargout{4} = z_vec;
varargout{5} = FINAL_CORRECTED_CENTROID;
varargout{6} = grad_out_LUMP_UPDATE;
varargout{7} = missing_idx;
varargout{8} = FIB_PIX_CORD;

end