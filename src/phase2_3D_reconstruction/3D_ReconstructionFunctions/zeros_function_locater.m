function [blocks] = zeros_function_locater(non_rogue_area_find)
% This function identifies contiguous non-zero blocks within a 1D array.
% It is typically used to locate continuous regions (e.g., fibers or valid data)
% within a binary or logical array indicating non-rogue areas.
%
% INPUT:
%   non_rogue_area_find : 1D array (logical or numeric) where non-zero indicates valid region
%
% OUTPUT:
%   blocks : Cell array, each containing indices of a contiguous non-zero block

% Multiply the logical/non-zero array by its linear indices to get actual positions
rogue_area_find_idx  = double(non_rogue_area_find) .* (1:numel(non_rogue_area_find));

% Pad with zeros at both ends to detect boundaries at edges
wrap = [0, rogue_area_find_idx, 0];

% Identify changes from zero to non-zero (start of block) and non-zero to zero (end of block)
temp = diff(wrap ~= 0);
blockStart = find(temp == 1) + 1;   % Indices where blocks start
blockEnd   = find(temp == -1);      % Indices where blocks end

% Extract contiguous non-zero blocks into a cell array
blocks = arrayfun(@(bId) wrap(blockStart(bId):blockEnd(bId)), ...
                  1:numel(blockStart), 'UniformOutput', false);

% Optional (commented out) for getting just first and last index of each block
% pos_neg_temp = blocks;
% chunks_fibers_idx_cell = cellfun(@(x)[x(1) x(end)], pos_neg_temp, 'uni', 0);

end