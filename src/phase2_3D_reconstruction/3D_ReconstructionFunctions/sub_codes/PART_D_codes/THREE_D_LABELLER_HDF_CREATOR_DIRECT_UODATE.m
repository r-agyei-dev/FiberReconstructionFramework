function [varargout] = THREE_D_LABELLER_HDF_CREATOR_DIRECT_UODATE(varargin)

%% ===================== Load Input Variables =====================
% Input arguments (in order):
%   1. check_mat_center_GLOBAL          - Cell array with linear indices of fibers
%   2. size_length                      - Dimensions of the 3D volume [X Y Z]
%   3. lump                             - Flag indicating whether to lump fibers together
%   4. hdfsave                          - Flag to save output in HDF/XDMF format
%   5. Centroid_MID_cell_LUMP_UPDATE    - Centroid coordinates of fibers
%   6. load_variable                    - Determines voxel ordering (1: XYZ, 2: YZX, 3: ZXY)
%   7. FINAL_CORRECTED_CENTROID         - Corrected centroid coordinates
%   8. idx                               - Index for labeling files
%   9. save_path                         - Directory to save HDF/XDMF files
%   10. sieve_num                        - Identifier for HDF/XDMF dataset

check_mat_center_GLOBAL            =  varargin{1};
size_length                        =  varargin{2};
lump                               =  varargin{3};
hdfsave                            =  varargin{4};
Centroid_MID_cell_LUMP_UPDATE      =  varargin{5};
load_variable                      =  varargin{6};
FINAL_CORRECTED_CENTROID           =  varargin{7};
idx                                =  varargin{8};
save_path                          =  varargin{9};
sieve_num                          =  varargin{10};

%% ===================== Remove Empty Fibers =====================
% Eliminate any empty entries in fiber indices or centroids
Centroid_MID_cell_LUMP_UPDATE(cell2mat(cellfun(@(x)isempty(x),check_mat_center_GLOBAL,'uni',0))) = [];
FINAL_CORRECTED_CENTROID(cell2mat(cellfun(@(x)isempty(x),check_mat_center_GLOBAL,'uni',0))) = [];
check_mat_center_GLOBAL(cell2mat(cellfun(@(x)isempty(x),check_mat_center_GLOBAL,'uni',0))) = [];

%% ===================== Prepare Indexing for Height =====================
if lump == 1
    % Combine all fibers into a single vector for height calculation
    check_mat_center_GLOBAL_VEC = [check_mat_center_GLOBAL{:}];
    Idx_Height = [];
else
    % Keep fibers separate for height indexing
    check_mat_center_GLOBAL_VEC = check_mat_center_GLOBAL;
    Idx_Height = cell2mat(cellfun(@(x)size(x,2),Centroid_MID_cell_LUMP_UPDATE,'uni',0))';
end

% Assign a second column for fiber numbering and flip columns for indexing
Idx_Height(:,2) = 1:size(Idx_Height(:,1),1);
Idx_Height = fliplr(Idx_Height);

%% ===================== Initialize 3D Volume and Fiber Arrays =====================
% TEMP_MAT stores the 3D labeled volume
% FIB_LIN_IDX stores linear indices for each fiber
% z_vec stores the range of z-coordinates for each fiber
TEMP_MAT = zeros([max(size_length) max(size_length) min(size_length)]);
FIB_LIN_IDX = cell(1,length(check_mat_center_GLOBAL_VEC));
z_vec = cell(1,length(check_mat_center_GLOBAL_VEC));

%% ===================== Populate 3D Volume =====================
for ii = 1:length(check_mat_center_GLOBAL_VEC)
    % Convert linear indices to subscript coordinates
    [x,y,z] = ind2sub(size_length, check_mat_center_GLOBAL_VEC{ii});
    z_vec{ii} = z(1):z(end);
    xyz = [x, y, z];
    
    % Determine permutation of axes based on load_variable
    if load_variable == 2       % YZX ordering
        perm_idx = [3 1 2];  
    elseif load_variable == 3   % ZXY ordering
        perm_idx = [2 3 1];
    end
    
    if load_variable > 1
        size_length_p = size_length(perm_idx);
        xyz = xyz(:, perm_idx);
    else
        size_length_p = size_length;   
    end
    
    % Convert to linear indices in permuted volume and fill TEMP_MAT
    FIB_LIN_IDX{ii} = sub2ind(size_length_p, xyz(:,1), xyz(:,2), xyz(:,3));
    TEMP_MAT(FIB_LIN_IDX{ii}) = ii;    
end

%% ===================== HDF/XDMF Saving =====================
if hdfsave == 1
    % Create filenames based on load_variable and sieve_num
    if load_variable == 1
        nom_bin  = sprintf('%s%d%s','REFINED_FIBERS_xyz_',sieve_num,'.h5');
        nom_file = sprintf('%s%d%s','REFINED_FIBERS_xyz_',sieve_num,'.xdmf');
    elseif load_variable == 2
        nom_bin  = sprintf('%s%d%s','REFINED_FIBERS_yzx_',sieve_num,'.h5');
        nom_file = sprintf('%s%d%s','REFINED_FIBERS_yzx_',sieve_num,'.xdmf');
    elseif load_variable == 3
        nom_bin  = sprintf('%s%d%s','REFINED_FIBERS_zxy_',sieve_num,'.h5');
        nom_file = sprintf('%s%d%s','REFINED_FIBERS_zxy_',sieve_num,'.xdmf');
    end
    
    % Delete existing files before writing
    delete(nom_bin);
    
    % Create HDF5 dataset with chunking and compression
    chunkSize = [200,200,50];
    h5create(nom_bin,'/CELL_DATA/FIBER_RECON',size(TEMP_MAT),'Datatype','double','ChunkSize',chunkSize,'Deflate',6);
    
    [nx,ny,nz] = size(TEMP_MAT);
    dz = chunkSize(3);
    zStart = 1:dz:nz;
    zEnd = dz:dz:nz;
    if zStart(end) == nz
        zStart = zStart(1:end-1);
    end
    if zEnd(end) ~= nz
        zEnd = [zEnd, nz];
    end
    
    % Write the data in z-chunks to HDF5
    for i = 1:numel(zStart)
        dz = zEnd(i) - zStart(i) + 1;
        if dz > 0
            start = [1 1 zStart(i)];        % Start of dataset
            count = [nx ny dz];             % Block size to write
            data = TEMP_MAT(:,:,zStart(i):zEnd(i)); % Extract data slice
            h5write(nom_bin,'/CELL_DATA/FIBER_RECON',data,start,count);
        else
            disp('dz == 0, cannot write data')
            break
        end
    end
    
    % Prepare attributes for XDMF
    h5att_name{1} = sprintf('%s%d','Reconstructed_Labels_check',idx);
    h5string{1}   = sprintf('%s%s',nom_bin,':/CELL_DATA/FIBER_RECON');
    
    % Generate XDMF file
    file_name = nom_file;
    z_depth = size(TEMP_MAT,3);
    x_breadth = size(TEMP_MAT,1);
    y_breadth = size(TEMP_MAT,2);
    
    XDMF_TEXT_GENERATOR(file_name,h5string,h5att_name,z_depth,y_breadth,x_breadth);
    
    % Move files to save_path
    movefile(nom_bin, save_path);
    movefile(nom_file, save_path);
end

%% ===================== Assign Outputs =====================
% 1. TEMP_MAT             - 3D labeled fiber volume
% 2. Idx_Height           - Fiber height and index information
% 3. FIB_LIN_IDX          - Linear indices for each fiber
% 4. z_vec                - Z-coordinate ranges for each fiber
% 5. FINAL_CORRECTED_CENTROID - Centroid coordinates of fibers

varargout{1} = TEMP_MAT;
varargout{end+1} = Idx_Height;
varargout{end+1} = FIB_LIN_IDX;
varargout{end+1} = z_vec;
varargout{end+1} = FINAL_CORRECTED_CENTROID;

end