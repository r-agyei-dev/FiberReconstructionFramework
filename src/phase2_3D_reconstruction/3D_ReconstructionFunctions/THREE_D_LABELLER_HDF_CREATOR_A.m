% This function generates a 3D labeled volume of fibers from linear indices.
% It can also compute the height indices, optionally save as HDF/XDMF, 
% and account for potential flipping or broken fibers.

function [TEMP_MAT,Idx_Height] = THREE_D_LABELLER_HDF_CREATOR_A(varargin)

% INPUT VARIABLES:
% Linear_Index_center_GLOBAL       = Cell array of linear indices representing fibers in 3D volume
% size_length                      = Size of the 3D reconstruction volume [rows, cols, slices]
% lump                             = Boolean flag indicating if fiber data is lumped together (1 = lumped)
% hdfsave                          = Boolean flag to create HDF/XDMF files
% Centroid_MID_cell_LUMP_UPDATE    = Cell array containing centroid pixel locations of fibers
% flip_query                       = Boolean flag to indicate if 3D volume is flipped
% load_variable                    = Axis orientation of serial sections (1=x,2=y,3=z)
% save_path                        = Directory path to save HDF/XDMF files
% new_idx                          = Index of the current tomogram dataset

Linear_Index_center_GLOBAL       =  varargin{1};
size_length                      =  varargin{2};
lump                             =  varargin{3}; 
hdfsave                          =  varargin{4};
Centroid_MID_cell_LUMP_UPDATE    =  varargin{5};
flip_query                       =  varargin{6};
load_variable                    =  varargin{7};
save_path                        =  varargin{8};
new_idx                          =  varargin{9};

% -------------------- PROCESS FIBER INDICES --------------------
if lump == 1
    % If fibers are lumped, combine all linear indices into a single vector
    check_mat_center_GLOBAL_VEC = [Linear_Index_center_GLOBAL{:}];
    Idx_Height = [];  % No height info when lumped
else 
    % If fibers are separate, retain original indices and compute height info
    check_mat_center_GLOBAL_VEC = Linear_Index_center_GLOBAL;
    Idx_Height = cellfun(@(x)size(x,1),Centroid_MID_cell_LUMP_UPDATE,'uni',0);
end

% Initialize a 3D matrix to store fiber labels
TEMP_MAT = zeros(size_length);

% Loop through each fiber and label its voxels in TEMP_MAT
for ii = 1:length(check_mat_center_GLOBAL_VEC)
   TEMP_MAT(check_mat_center_GLOBAL_VEC{ii}) = ii;    
end

%% ----------------------------------------------------------------
%% OPTIONAL VISUALIZATION AND FILE SAVING (commented out)
% TEMP_MAT_VISUAL = TEMP_MAT; 
% 
% % If 3D volume is flipped, permute axes according to load_variable
% if flip_query == 1
%      if load_variable == 3  % yzx orientation
%          TEMP_MAT_VISUAL = permute(TEMP_MAT_VISUAL,[3 1 2]);   
%      elseif load_variable == 2  % zxy orientation
%          TEMP_MAT_VISUAL = permute(TEMP_MAT_VISUAL,[2 3 1]);
%      end
% end
%
% % If saving is enabled, create HDF and XDMF files for visualization
% if hdfsave == 1
%     
%     % Depending on axis orientation, assign filenames and save volume
%     if load_variable == 1  % xyz orientation
%         nom_bin = sprintf('RECON_xyz_labeller_%d_.h5', new_idx);
%         file_name = sprintf('RECON_xyz_labeller_%d_.xdmf', new_idx);
%         delete(nom_bin)
%         h5create(nom_bin,'/CELL_DATA/FIBER_RECON',size(TEMP_MAT_VISUAL),'Datatype','double');
%         h5write(nom_bin,'/CELL_DATA/FIBER_RECON',double(TEMP_MAT_VISUAL));
%         h5att_name{1} = 'Labels_xyz';
%         h5string{1}   = sprintf('%s:/CELL_DATA/FIBER_RECON', nom_bin);
%     elseif load_variable == 2  % zxy orientation
%         % Similar saving logic for zxy
%     elseif load_variable == 3  % yzx orientation
%         % Similar saving logic for yzx
%     end
%     
%     % Generate XDMF metadata for visualization
%     z_depth = size(TEMP_MAT_VISUAL,3);
%     x_breadth = size(TEMP_MAT_VISUAL,1);
%     y_breadth = size(TEMP_MAT_VISUAL,2);
%     XDMF_TEXT_GENERATOR(file_name,h5string,h5att_name,z_depth,y_breadth,x_breadth)
% 
%     % Move generated files to save_path
%     movefile(nom_bin,save_path)
%     movefile(file_name,save_path)
% end