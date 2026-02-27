function [varargout] = THREE_D_LABELLER_HDF_CREATOR_DIRECT_UODATE_v2_TT(varargin)

% =========================================================================
% THREE_D_LABELLER_HDF_CREATOR_DIRECT_UODATE_v2_TT
%
% PURPOSE:
%   Processes refined fiber data to:
%     1) Clean and validate fiber index/centroid data
%     2) Correct split fibers
%     3) Recompute fiber orientation angles
%     4) Reconstruct a labeled 3D volume
%     5) Export the labeled volume to HDF format
%
% INPUTS (via varargin):
%   1  - Linear_Index_Final           : Cell array of fiber linear indices
%   2  - orig_size_length             : Original volume size [X Y Z]
%   3  - size_length                  : Working volume size
%   4  - hdfsave                      : Flag for HDF saving
%   5  - Pixel_Idx_Per_Slice          : Pixel indices per slice
%   6  - load_variable                : Axis ordering flag
%   7  - FINAL_CORRECTED_CENTROID     : Fiber centroid coordinates
%   8  - tomo_dataset_idx             : Dataset index for HDF
%   9  - save_path                    : Output directory
%   10 - sieve_indicator              : Flag for sieved vs unsieved fibers
%   11 - angular_thresh               : (reserved/unused here)
%   12 - height_variable_MAIN         : Fiber height info
%   13 - crop_region                  : Z cropping bounds
%
% OUTPUTS:
%   1  - TEMP_MAT                     : 3D labeled fiber volume
%   2  - Idx_Height                   : Fiber height table
%   3  - FIB_LIN_IDX                  : Linear indices per fiber
%   4  - z_vec                        : Z-span per fiber
%   5  - FINAL_CORRECTED_CENTROID_LUMP
%   6  - grad_out_LUMP_UPDATE         : Out-of-plane angles
%   7  - missing_idx                  : Missing fiber indices
%   8  - FIB_PIX_CORD                 : Per-fiber voxel coordinates
% =========================================================================

Linear_Index_Final                 =  varargin{1};
orig_size_length                   =  varargin{2};
size_length                        =  varargin{3};
hdfsave                            =  varargin{4};
Pixel_Idx_Per_Slice                =  varargin{5};
load_variable                      =  varargin{6};
FINAL_CORRECTED_CENTROID           =  varargin{7};
tomo_dataset_idx                   =  varargin{8};
save_path                          =  varargin{9};
sieve_indicator                    =  varargin{10};
angular_thresh                     =  varargin{11};
height_variable_MAIN               =  varargin{12};
crop_region                        =  varargin{13};

% -------------------------------------------------------------------------
% Select output naming convention based on sieving flag
% -------------------------------------------------------------------------
if sieve_indicator == 1
    % Use naming for sieved (quality-filtered) fibers
    nom_bin_prefix            =   'REFINED_FIBERS_SIEVED';
    nom_file_prefix           =   'REFINED_FIBERS_SIEVED';
else
    % Use naming for unsieved fibers
    nom_bin_prefix            =   'REFINED_FIBERS_UNSIEVED';
    nom_file_prefix           =   'REFINED_FIBERS_UNSIEVED';    
end

% -------------------------------------------------------------------------
% Remove empty entries across all synchronized fiber containers
% -------------------------------------------------------------------------
Pixel_Idx_Per_Slice(cell2mat(cellfun(@(x)isempty(x),Linear_Index_Final,'uni',0)))                =   [];
FINAL_CORRECTED_CENTROID((cell2mat(cellfun(@(x)isempty(x),Linear_Index_Final,'uni',0))))         =   [];
height_variable_MAIN((cell2mat(cellfun(@(x)isempty(x),Linear_Index_Final,'uni',0))))             =   [];
Linear_Index_Final((cell2mat(cellfun(@(x)isempty(x),Linear_Index_Final,'uni',0))))               =   [];

%% ------------------------------------------------------------------------
% Investigate and correct fiber splits in the refined microstructure
% ------------------------------------------------------------------------
disp('Investigate the splits within the refined microstructure')

minor_axis_cell    =   [];

% Split correction and relumping
[Linear_Index_Final_LUMP,~,FINAL_CORRECTED_CENTROID_LUMP,~,~] = ...
    SPLITTER_CORRECTOR_FUNC_CODE(Linear_Index_Final,minor_axis_cell,...
    FINAL_CORRECTED_CENTROID,height_variable_MAIN,Pixel_Idx_Per_Slice);

FINAL_CORRECTED_CENTROID_LUMP = Splitter_Lumper_func(FINAL_CORRECTED_CENTROID_LUMP);                                                           
Linear_Index_Final_LUMP       = Splitter_Lumper_func(Linear_Index_Final_LUMP);    

%% ------------------------------------------------------------------------
% Remove any empty fibers after lumping
% ------------------------------------------------------------------------
disp('delete the empty cells in the check_mat_center_GLOBAL_VEC or the FINAL_CORRECTED_CENTROID')

missing_idx_1 = cell2mat(cellfun(@(x)isempty(x),Linear_Index_Final_LUMP,'uni',0));
missing_idx_2 = cell2mat(cellfun(@(x)isempty(x),FINAL_CORRECTED_CENTROID_LUMP,'uni',0));

Linear_Index_Final_LUMP(missing_idx_1 | missing_idx_2)     =  [];
FINAL_CORRECTED_CENTROID_LUMP(missing_idx_1 | missing_idx_2) =  [];

%% ------------------------------------------------------------------------
% Compute fiber orientation angles
% ------------------------------------------------------------------------
disp('computing the angles')

% Identify single-point fibers
single_temp = cell2mat(cellfun(@(x)size(x,1),...
                FINAL_CORRECTED_CENTROID_LUMP,'uni',0)) == 1;

% In-plane angle (row/column plane)
in_plane_cell = cellfun(@(x) ...
    atand((x(1,2)-x(end,2))/(x(end,1)-x(1,1))),...
    FINAL_CORRECTED_CENTROID_LUMP,'uni',0);

% Force single-point fibers to zero angle
[in_plane_cell{single_temp}] = deal(0);

disp('For multiple fibers protruding slices (make zero all the NaN values of the in-plane_cell angles') 
tic

% Handle zero-length projections
[in_plane_cell{cell2mat(cellfun(@(x) ...
    norm(x(1,1:2)-x(end,1:2)),...
    FINAL_CORRECTED_CENTROID_LUMP,'uni',0)) == 0}] = deal(0);

% Out-of-plane angle calculation
out_plane_angle = cellfun(@(x,p) ...
    90 - acosd((x(end,1)-x(1,1)) / ...
    (norm(x(end,:)-x(1,:))*cosd(p))),...
    FINAL_CORRECTED_CENTROID_LUMP,in_plane_cell,'uni',0);
toc

% Replace NaNs with zero
disp('for all the NaN make them 90 degrees as well')
[out_plane_angle{cell2mat(cellfun(@(x)isnan(x),out_plane_angle,'uni',0))}] = deal(0);
out_plane_angle = cell2mat(out_plane_angle);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Prepare height and angle lookup tables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('UPDATE THE VARIABLES FOR THE OUT-OF-PLANE ANGLE AND THE IDX_HEIGHT')
tic
grad_out_LUMP_UPDATE   = out_plane_angle;

Idx_Height             = cell2mat(cellfun(@(x)size(x,1),...
                                FINAL_CORRECTED_CENTROID_LUMP,'uni',0))';
Idx_Height(:,2)        = 1:size(Idx_Height(:,1),1);
Idx_Height             = fliplr(Idx_Height);
toc

%% ------------------------------------------------------------------------
% Reconstruct labeled 3D volume
% ------------------------------------------------------------------------
disp('creating the 3D volume')

crop_height  = orig_size_length(3);
TEMP_MAT     = zeros([orig_size_length]); % Memory-efficient allocation

FIB_LIN_IDX  = cell(1,length(Linear_Index_Final_LUMP));
z_vec        = cell(1,length(Linear_Index_Final_LUMP));
FIB_PIX_CORD = cell(1,length(Linear_Index_Final_LUMP));

tic
for ii = 1:length(Linear_Index_Final_LUMP)

    % Convert linear indices to subscripts
    [x,y,z] = ind2sub(size_length,Linear_Index_Final_LUMP{ii});
    z_vec{ii} = z(1):z(end);
    xyz = [x,y,z];
    
    % -------------------------------------------------------------
    % Handle different data loading permutations
    % -------------------------------------------------------------
    if load_variable == 2          % zxy ordering
        perm_idx = [2 3 1];  
    elseif load_variable == 3      % yzx ordering
        perm_idx = [3 1 2];
    end     
     
    if load_variable > 1
        size_length_p = size_length(perm_idx);
        xyz = xyz(:,perm_idx);
    else 
        size_length_p = size_length;   
    end
     
    size_length_p = [size_length_p(1:2) crop_height];
     
    % -------------------------------------------------------------
    % Apply crop height normalization
    % -------------------------------------------------------------
    height_leveller = crop_region(1) + 1;
    xyz = [xyz(:,1:2)  xyz(:,3) - height_leveller*ones(size(xyz,1),1)];
      
    % Remove voxels outside cropped bounds
    xyz(xyz(:,3)<1,:) = [];
    xyz(xyz(:,3)>crop_height,:) = [];
      
    % Convert back to linear indices and label volume
    FIB_LIN_IDX{ii} = sub2ind(size_length_p,xyz(:,1),xyz(:,2),xyz(:,3));
    TEMP_MAT(FIB_LIN_IDX{ii}) = ii;  
      
    % Store voxel coordinates per fiber
    FIB_PIX_CORD{ii} = xyz;
end
toc

% -------------------------------------------------------------------------
% (Reserved) Missing fiber detection
% -------------------------------------------------------------------------
missing_idx = [];

%% ------------------------------------------------------------------------
% Export labeled volume to HDF
% ------------------------------------------------------------------------
disp('creating hdf files')
tic
[~] = MAIN_HDF_CREATOR(TEMP_MAT,nom_bin_prefix,nom_file_prefix,...
        tomo_dataset_idx,hdfsave,load_variable,save_path);
toc

%% ------------------------------------------------------------------------
% Outputs
% ------------------------------------------------------------------------
varargout{1}      = TEMP_MAT;
varargout{end+1}  = Idx_Height;
varargout{end+1}  = FIB_LIN_IDX;
varargout{end+1}  = z_vec;
varargout{end+1}  = FINAL_CORRECTED_CENTROID_LUMP;
varargout{end+1}  = grad_out_LUMP_UPDATE;
varargout{end+1}  = missing_idx;
varargout{end+1}  = FIB_PIX_CORD;

end