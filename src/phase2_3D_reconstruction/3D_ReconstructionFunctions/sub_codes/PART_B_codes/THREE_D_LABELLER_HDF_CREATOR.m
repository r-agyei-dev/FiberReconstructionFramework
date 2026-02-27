% This function generates 
% (a) the volume image for the respective planes  
% (b) the Idx_Height for the respective planes 
% (c) writes an hdf file if the hdfsave variable is invoked 
% (d) flips the image volume if reconstruction plane is not along the z axis 
% (e) update has been created to fix the potential broken fiber fix in the code as well

%% Version Update: 09/15/2018

function [TEMP_MAT,Idx_Height] = THREE_D_LABELLER_HDF_CREATOR(varargin)
Linear_Index_center_GLOBAL       =  varargin{1};
size_length                      =  varargin{2};
lump                             =  varargin{3}; 
hdfsave                          =  varargin{4};
Centroid_MID_cell_LUMP_UPDATE    =  varargin{5};
flip_query                       =  varargin{6};
load_variable                    =  varargin{7};
save_path                        =  varargin{8};
new_idx                          =  varargin{9};
fiber_idx_split_amt_VARIABLE     =  varargin{10};
% TEMP_CRUDE                       =  varargin{11};

if lump == 1
check_mat_center_GLOBAL_VEC = [Linear_Index_center_GLOBAL{:}];
Idx_Height = [];
else 
check_mat_center_GLOBAL_VEC = Linear_Index_center_GLOBAL;
Idx_Height = cellfun(@(x)size(x,1),Centroid_MID_cell_LUMP_UPDATE,'uni',0);
end

TEMP_MAT = zeros(size_length);
for ii = 1:length(check_mat_center_GLOBAL_VEC)
   TEMP_MAT(check_mat_center_GLOBAL_VEC{ii}) = ii;    
end

% TEMP_CRUDE(TEMP_MAT~=0) = 0;
%%
TEMP_MAT_VISUAL = TEMP_MAT; 

if flip_query == 1
     if load_variable == 3                        % refers to yzx
         TEMP_MAT_VISUAL = permute(TEMP_MAT_VISUAL,[3 1 2]);   
     elseif load_variable == 2                    % refers to zxy
         TEMP_MAT_VISUAL = permute(TEMP_MAT_VISUAL,[2 3 1]);
     end
end
%%

if hdfsave == 1
    
    if load_variable == 1

nom_bin = sprintf('%s%d%s','RECON_xyz_labeller_',new_idx,'_.h5');
file_name = sprintf('%s%d%s','RECON_xyz_labeller_',new_idx,'_.xdmf');
delete(nom_bin)
h5create(sprintf('%s',nom_bin),'/CELL_DATA/FIBER_RECON',size(TEMP_MAT_VISUAL),'Datatype','double');
h5write(sprintf('%s',nom_bin),'/CELL_DATA/FIBER_RECON',double(TEMP_MAT_VISUAL));
h5att_name{1}  =  sprintf('%s','Labels_xyz');
h5string{1}    =  sprintf('%s%s',nom_bin,':/CELL_DATA/FIBER_RECON');


% fiber_idx_split_amt_VARIABLE_xyz = fiber_idx_split_amt_VARIABLE;
% save(fullfile(save_path,'fiber_idx_split_amt_VARIABLE_xyz.mat'),'fiber_idx_split_amt_VARIABLE_xyz','-v7.3');

       
    elseif load_variable == 2

nom_bin = sprintf('%s%d%s','RECON_zxy_labeller_',new_idx,'_.h5');
file_name = sprintf('%s%d%s','RECON_zxy_labeller_',new_idx,'_.xdmf');
delete(nom_bin)
h5create(sprintf('%s',nom_bin),'/CELL_DATA/FIBER_RECON',size(TEMP_MAT_VISUAL),'Datatype','double');
h5write(sprintf('%s',nom_bin),'/CELL_DATA/FIBER_RECON',double(TEMP_MAT_VISUAL));
h5att_name{1}  =  sprintf('%s','Labels_zxy');
h5string{1}    =  sprintf('%s%s',nom_bin,':/CELL_DATA/FIBER_RECON');


% fiber_idx_split_amt_VARIABLE_zxy = fiber_idx_split_amt_VARIABLE;
% save(fullfile(save_path,'fiber_idx_split_amt_VARIABLE_zxy.mat'),'fiber_idx_split_amt_VARIABLE_zxy','-v7.3');

   elseif load_variable == 3
        
nom_bin = sprintf('%s%d%s','RECON_yzx_labeller_',new_idx,'_.h5');
file_name = sprintf('%s%d%s','RECON_yzx_labeller_',new_idx,'_.xdmf');
delete(nom_bin)
h5create(sprintf('%s',nom_bin),'/CELL_DATA/FIBER_RECON',size(TEMP_MAT_VISUAL),'Datatype','double');
h5write(sprintf('%s',nom_bin),'/CELL_DATA/FIBER_RECON',double(TEMP_MAT_VISUAL));
h5att_name{1}  =  sprintf('%s','Labels_yzx');
h5string{1}    =  sprintf('%s%s',nom_bin,':/CELL_DATA/FIBER_RECON');


% fiber_idx_split_amt_VARIABLE_yzx = fiber_idx_split_amt_VARIABLE;
% save(fullfile(save_path,'fiber_idx_split_amt_VARIABLE_yzx.mat'),'fiber_idx_split_amt_VARIABLE_yzx','-v7.3');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end

% % ========================================================================

z_depth = size(TEMP_MAT_VISUAL,3);
x_breadth = size(TEMP_MAT_VISUAL,1);
y_breadth = size(TEMP_MAT_VISUAL,2);

XDMF_TEXT_GENERATOR(file_name,h5string,h5att_name,z_depth,y_breadth,x_breadth)

movefile(nom_bin,save_path)
movefile(file_name,save_path)

end


