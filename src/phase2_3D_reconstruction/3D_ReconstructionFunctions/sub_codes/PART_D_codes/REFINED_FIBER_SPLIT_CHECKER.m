% REFINED_FIBER_SPLIT_CHECKER
% This function identifies and resolves splits within fiber blocks.
% It operates on a cell array of fiber centroid positions and their corresponding
% linear indices, separating fibers into contiguous blocks based on shared pixels.
% 
% Inputs:
%   Centroid_MID_cell_LUMP_UPDATE    - Cell array of fiber centroid positions [Nx3]
%   check_mat_center_GLOBAL_VEC      - Cell array of linear indices corresponding to fiber centroids
%   FINAL_CORRECTED_CENTROID         - Cell array of corrected fiber centroid positions
%   size_length                      - Size of the 3D volume [rows, columns, slices]
%
% Outputs:
%   varargout{1} - Concatenated linear indices after splitting fiber blocks
%   varargout{2} - Concatenated fiber centroids after splitting fiber blocks

function [varargout] = REFINED_FIBER_SPLIT_CHECKER(varargin)

Centroid_MID_cell_LUMP_UPDATE    =   varargin{1};
check_mat_center_GLOBAL_VEC      =   varargin{2};
FINAL_CORRECTED_CENTROID         =   varargin{3};
size_length                      =   varargin{4};

%% STEP 1: Identify splits within fiber blocks
% Initialize a cell array to store contiguous blocks for each fiber
blocks_cell = cell(1, size(Centroid_MID_cell_LUMP_UPDATE, 2));

for kk = 1:size(Centroid_MID_cell_LUMP_UPDATE,2)
    
    % Convert each fiber's centroid coordinates into linear indices
    aa = cellfun(@(x)sub2ind(size_length(1:2), x(:,2), x(:,1)), Centroid_MID_cell_LUMP_UPDATE{kk}, 'uni', 0);
    num_el = 1:size(aa, 2);
    
    % Detect overlaps between consecutive fibers
    split_indicator = zeros(1, length(num_el)-1);
    for ii = 1:length(num_el)-1
        split_indicator(ii) = nnz(ismember(aa{ii}, aa{ii+1}));
    end
    
    % Pad with zeros to simplify detection of block starts and ends
    wrap = [0, [double(split_indicator ~= 0) 1] .* (1:length(num_el)), 0];
    temp = diff(wrap ~= 0);
    blockStart = find(temp == 1) + 1;
    blockEnd   = find(temp == -1);
    
    % Extract contiguous blocks of fibers
    blocks = arrayfun(@(bId) wrap(blockStart(bId):blockEnd(bId)), 1:numel(blockStart), 'UniformOutput', false);
    
    % Append the end of each block if necessary
    for tt = 1:size(blocks,1)
        if blocks{tt} ~= length(num_el)
            blocks{tt}(end+1) = blocks{tt}(end)+1;
        end
    end
    
    % Include isolated fragments (fibers not part of any block)
    blocks_cell{kk} = horzcat(blocks, num2cell(num_el(~ismember(num_el, horzcat(blocks{:}))),1));
end   

%% STEP 2: Split linear indices and centroids based on detected blocks
% Identify fiber blocks that contain multiple fragments
split_idx = find(cell2mat(cellfun(@(x)size(x,2) > 1, blocks_cell, 'uni', 0)));

for bb = 1:numel(split_idx)
    
    % Assign each pixel index a unique multiplier to track fiber membership
    premultiplier_cell = num2cell(1:size(Centroid_MID_cell_LUMP_UPDATE{split_idx(bb)}, 2),1);
    ones_vector_cell = cellfun(@(t)ones(1,size(t,1)), Centroid_MID_cell_LUMP_UPDATE{split_idx(bb)}, 'uni',0);
    num_tmp = cellfun(@(x,y)x*y, premultiplier_cell, ones_vector_cell, 'uni', 0);
    num_tmp = horzcat(num_tmp{:})';
    
    % Split linear indices according to the detected blocks
    check_mat_tmp = cellfun(@(x)check_mat_center_GLOBAL_VEC{split_idx(bb)}(ismember(num_tmp, x)), blocks_cell{split_idx(bb)}, 'uni', 0);
    check_mat_center_GLOBAL_VEC{split_idx(bb)} = check_mat_tmp;
    
    % Split centroid positions according to the detected blocks
    fiber_centroid_tmp = cellfun(@(x)FINAL_CORRECTED_CENTROID{split_idx(bb)}(ismember(FINAL_CORRECTED_CENTROID{split_idx(bb)}(:,3), x), :), blocks_cell{split_idx(bb)}, 'uni',0);
    FINAL_CORRECTED_CENTROID{split_idx(bb)} = fiber_centroid_tmp;
end

%% STEP 3: Concatenate results for output
varargout{1} = [check_mat_center_GLOBAL_VEC{:}];
varargout{2} = [FINAL_CORRECTED_CENTROID{:}];

end