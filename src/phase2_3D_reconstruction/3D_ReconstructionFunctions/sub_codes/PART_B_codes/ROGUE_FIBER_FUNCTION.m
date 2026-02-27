function [varargout] = ROGUE_FIBER_FUNCTION(varargin)

% ROGUE_FIBER_FUNCTION
% -------------------------------------------------------------------------
% Identifies and groups "rogue" fiber segments, reconstructs their voxel
% representation in 3D, and exports the labeled volume to HDF5/XDMF for
% visualization.
%
% INPUTS (via varargin)
%   1  LUMP_ROGUE            : cell array of rogue slice indices per fiber
%   2  ROGUE_BEGIN_SLICE     : starting slice index for each rogue fiber
%   3  idx_after_first_thresh: original fiber indices
%   4  main_cell             : main fiber structure containing Pixel_List
%   5  save_path             : destination folder for output files
%   6  size_length           : size of the 3D volume
%   7  new_idx               : identifier used in output filenames
%   8  flip_query            : flag to permute volume for visualization
%   9  load_variable         : orientation selector for permutation
%
% OUTPUTS
%   Structure_Rogue : structure describing rogue chunks
%   TEMP_ROGUE      : labeled rogue volume
%   rogue_ends      : start/end slice indices of each rogue chunk
% -------------------------------------------------------------------------

LUMP_ROGUE               =   varargin{1};
ROGUE_BEGIN_SLICE        =   varargin{2};
idx_after_first_thresh   =   varargin{3};
main_cell                =   varargin{4};
save_path                =   varargin{5};
size_length              =   varargin{6};
new_idx                  =   varargin{7};
flip_query               =   varargin{8};
load_variable            =   varargin{9};

% Preallocate containers for separated rogue chunks and metadata
chunks_fibers_idx_cell            =  cell(1,size(LUMP_ROGUE,2)); % chunked rogue indices
chunks_fibers_idx_cell_numel      =  cell(1,size(LUMP_ROGUE,2)); % length of each chunk
chunks_fibers_idx_cell_slice      =  cell(1,size(LUMP_ROGUE,2)); % starting slice per chunk
chunks_fibers_idx_cell_original   =  cell(1,size(LUMP_ROGUE,2)); % original fiber id

% -------------------------------------------------------------------------
% Split rogue indices into contiguous chunks for each fiber
% -------------------------------------------------------------------------
for pp = 1:size(LUMP_ROGUE,2)

if ~isempty(LUMP_ROGUE{pp})

% Create continuous span and logical mask of rogue positions
tmp_var_span = min(LUMP_ROGUE{pp}):max(LUMP_ROGUE{pp});
logical_var  = ismember(tmp_var_span,LUMP_ROGUE{pp}); 

% Label rogue positions
rogue_area_find_idx  =   double(logical_var).*(1:numel(logical_var));

% Pad with zeros to detect boundaries
wrap                 =   [0, rogue_area_find_idx, 0];

% Detect contiguous nonzero blocks
temp                 =   diff( wrap ~= 0 ) ;
blockStart           =   find( temp == 1 ) + 1 ;
blockEnd             =   find( temp == -1 ) ;
blocks               =   arrayfun(@(bId) wrap(blockStart(bId):blockEnd(bId)), ...
                                  1:numel(blockStart), 'UniformOutput', false);

% Convert block indices back to slice numbers
chunks_fibers_idx_cell{pp} = cellfun(@(x)tmp_var_span(x),blocks,'uni',0);

% Store chunk sizes
chunks_fibers_idx_cell_numel{pp} = ...
    cellfun(@(x)numel(x),chunks_fibers_idx_cell{pp},'uni',0);

% Record slice and original fiber metadata
chunks_fibers_idx_cell_slice{pp}    = ...
    ROGUE_BEGIN_SLICE(pp)*ones(1,length(chunks_fibers_idx_cell_numel{pp}));

chunks_fibers_idx_cell_original{pp} = ...
    idx_after_first_thresh(pp)*ones(1,length(chunks_fibers_idx_cell_numel{pp}));
end
                                                      
end 

% -------------------------------------------------------------------------
% Assemble rogue chunk structure
% -------------------------------------------------------------------------
Structure_Rogue.chunks_idx_LUMP          =   [chunks_fibers_idx_cell{:}];
Structure_Rogue.chunks_numel_LUMP        =   cell2mat([chunks_fibers_idx_cell_numel{:}]);
Structure_Rogue.chunk_begin_slice_LUMP   =   [chunks_fibers_idx_cell_slice{:}];
Structure_Rogue.chunks_original_LUMP     =   [chunks_fibers_idx_cell_original{:}];

%% ------------------------------------------------------------------------
% Build labeled 3D rogue volume
% ------------------------------------------------------------------------
TEMP_ROGUE = zeros(size_length);

for kk = 1:size(Structure_Rogue.chunks_idx_LUMP,2)
   if ~isempty(Structure_Rogue.chunks_idx_LUMP{kk})

       % Extract pixel lists for this rogue chunk
       PIXEL_extract = {main_cell{Structure_Rogue.chunks_original_LUMP(kk)} ...
                        (Structure_Rogue.chunks_idx_LUMP{kk}).Pixel_List};

       % Convert to linear indices in 3D volume
       IDX_extract = cellfun(@(x,y)sub2ind(size_length,x(:,2),x(:,1), ...
                     (y+Structure_Rogue.chunk_begin_slice_LUMP(kk)-1) ...
                     *ones(size(x(:,2),1),1)), ...
                     PIXEL_extract, ...
                     num2cell(Structure_Rogue.chunks_idx_LUMP{kk},1),'uni',0);
       
       % Determine z-extent of rogue segment
       [~,~,d] = ind2sub(size_length,vertcat(IDX_extract{:}));
       d = unique(d);
       rogue_ends{kk} = [d(1) d(end)];
       
       % Store voxel indices and label in volume
       Pix_IDX_vec{kk} = vertcat(IDX_extract{:});
       TEMP_ROGUE(Pix_IDX_vec{kk}) = -kk;
   end
clearvars PIXEL_extract   
end

Structure_Rogue.rogue_ends = rogue_ends;

% -------------------------------------------------------------------------
% Prepare visualization volume (optional permutation)
% -------------------------------------------------------------------------
TEMP_ROGUE_VISUAL = TEMP_ROGUE; 

if flip_query == 1
     if load_variable == 3            % orientation: yzx
         TEMP_ROGUE_VISUAL = permute(TEMP_ROGUE_VISUAL,[3 1 2]);   
     elseif load_variable == 2        % orientation: zxy
         TEMP_ROGUE_VISUAL = permute(TEMP_ROGUE_VISUAL,[2 3 1]);
     end
end

% -------------------------------------------------------------------------
% Export rogue volume to HDF5 + XDMF
% -------------------------------------------------------------------------
nom_bin  = sprintf('%s%d%s','ROGUE_vol',new_idx,'_.h5');
nom_file = sprintf('%s%d%s','ROGUE_vol',new_idx,'_.xdmf');

delete(nom_bin)

h5create(nom_bin,'/CELL_DATA/FIBER_RECON',size(TEMP_ROGUE_VISUAL),'Datatype','double');
h5write(nom_bin,'/CELL_DATA/FIBER_RECON',double(TEMP_ROGUE_VISUAL));

h5att_name{1}  =  'Rogue_Labels';
h5string{1}    =  sprintf('%s%s',nom_bin,':/CELL_DATA/FIBER_RECON');

z_depth   = size(TEMP_ROGUE_VISUAL,3);
x_breadth = size(TEMP_ROGUE_VISUAL,1);
y_breadth = size(TEMP_ROGUE_VISUAL,2);

XDMF_TEXT_GENERATOR(nom_file,h5string,h5att_name,z_depth,y_breadth,x_breadth)

movefile(nom_bin,save_path)
movefile(nom_file,save_path)

%% ------------------------------------------------------------------------
% Outputs
% ------------------------------------------------------------------------
varargout{1}       =   Structure_Rogue; 
varargout{end+1}   =   TEMP_ROGUE;
varargout{end+1}   =   rogue_ends;