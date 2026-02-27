% CELLFUN_CHECK
% This function analyzes a cell array of fiber centroid positions to identify
% blocks of consecutive fibers that share common pixels in their projected 2D slices.
% It returns the linear indices of each centroid and a cell array indicating 
% groups of fibers that are connected or separated.
%
% Inputs:
%   Centroid_MID_cell_LUMP_UPDATE - Cell array of updated fiber centroids [Nx3]
%   size_length                    - Size of the 3D volume [rows, columns, slices]
%
% Outputs:
%   aa        - Cell array containing linear indices for each fiber's 2D projection
%   sep_vec   - Cell array indicating contiguous blocks of fibers (connected by shared pixels)

function [aa, sep_vec] = CELLFUN_CHECK(Centroid_MID_cell_LUMP_UPDATE, size_length)

% Convert each centroid's 2D coordinates into linear indices for the 2D slice
aa = cellfun(@(x)sub2ind(size_length(1:2), x(:,2), x(:,1)), Centroid_MID_cell_LUMP_UPDATE, 'uni', 0);

% Initialize vector to record overlaps between consecutive fibers
split_indicator = zeros(1, size(aa, 2)-1);

% Check for shared pixels between consecutive fibers
for ii = 1 : size(aa, 2)-1
    split_indicator(ii) = nnz(ismember(aa{ii}, aa{ii+1}));
end

% If any overlap exists between fibers, identify contiguous blocks
if ~isempty(split_indicator ~= 0)
    
    % Pad the split indicator with zeros at the start and end to simplify block detection
    wrap = [0, [double(split_indicator ~= 0) 1] .* (1:size(aa,2)), 0];
    
    % Compute differences to find block starts (1) and block ends (-1)
    temp       = diff(wrap ~= 0);
    blockStart = find(temp == 1) + 1;
    blockEnd   = find(temp == -1);
    
    % Collect indices of each contiguous block of fibers
    blocks = arrayfun(@(bId) wrap(blockStart(bId):blockEnd(bId)), 1:numel(blockStart), 'UniformOutput', false);
    
    % Identify fibers that are completely isolated (no shared pixels)
    cell_tmp = num2cell(find(split_indicator == 0), 1);
    cell_tmp{end+1} = [];
    
    % Combine contiguous blocks with isolated fibers to form final separation vector
    sep_vec = cellfun(@(x, y)[x y], blocks, cell_tmp, 'uni', 0);
else
    sep_vec = [];
end

end