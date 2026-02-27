%% =============================================================================
% EPAD_FUNCTION_CODE_Fiber_Stitch_ver
% =============================================================================
% Description:
% This function identifies and processes overlapping 3D fiber indices across
% consecutive slices. It groups fibers that are contiguous in slice space,
% handles gaps, and prepares arrays for subsequent stitching or separation.
%
% INPUTS:
% -------
% curr_values : Cell array of current slice fiber indices
% succ_values : Cell array of successive slice fiber indices
% exponent    : Large multiplier to encode unique fiber pairs
%
% OUTPUTS:
% --------
% curr_succ_multiple_cell : Cell array of [current, successive] fiber matches per fiber
% split_array_consec_int  : Cell array indicating contiguous slice segments
% split_array_slice       : Cell array of slice indices for each contiguous segment
%
% NOTES:
% ------
% - Uses quotient/remainder trick to encode fiber matches uniquely
% - Handles non-contiguous slices and splits them into separate groups
% - Can be parallelized (parfor) for large datasets
%% =============================================================================

function [curr_succ_multiple_cell, split_array_consec_int, split_array_slice] = EPAD_FUNCTION_CODE_Fiber_Stitch_ver(curr_values, succ_values, exponent)

tic

%% ============================================================================
% STEP 1: Encode fiber pairs uniquely
% ============================================================================
% For each fiber in the current slice, combine with its successor using
% a large multiplier (exponent) to create unique identifiers
curr_succ_SUM = cell(1, numel(curr_values));
for ii = 1:numel(curr_values)
    curr_succ_SUM{ii} = unique(nonzeros(curr_values{ii}*exponent + succ_values{ii}));
end

% Concatenate all encoded values
UNIQUE_SUM = vertcat(curr_succ_SUM{:});

% Remove any values smaller than exponent (invalid or empty matches)
UNIQUE_SUM(UNIQUE_SUM < exponent) = [];

% Clear intermediate variables to save memory
disp('clear to relieve memory')
clearvars curr_succ_SUM

%% ============================================================================
% STEP 2: Decode fiber pairs
% ============================================================================
% Use integer division and remainder to extract original fiber indices
for ii = 1:numel(UNIQUE_SUM)
    quotient(ii)  = floor(UNIQUE_SUM(ii)/exponent); % current slice fiber
    remainder(ii) = mod(UNIQUE_SUM(ii), exponent); % successive slice fiber
end

curr_succ = [quotient' remainder'];
clearvars quotient remainder
toc

%% ============================================================================
% STEP 3: Initialize output containers
% ============================================================================
% Extract unique fibers from current slice
unique_fib_idx = nonzeros(unique(curr_succ(:,1)));

% Preallocate arrays for storing results
curr_succ_multiple_cell = cell(1, numel(unique_fib_idx));
split_array_consec_int  = cell(1, numel(unique_fib_idx));
split_array_slice       = cell(1, numel(unique_fib_idx));

%% ============================================================================
% STEP 4: Process each fiber individually
% ============================================================================
% For each unique fiber, extract all its matches and identify contiguous slices
tic
% parfor ii = 1:numel(unique_fib_idx) % Uncomment if more RAM is available
for ii = 1:numel(unique_fib_idx)
    
    % Extract matches for this fiber
    curr_succ_tmp = curr_succ;
    curr_succ_multiple = curr_succ_tmp(ismember(curr_succ_tmp(:,1), unique_fib_idx(ii)), :)';
    
    % Remove any zero slice entries
    curr_succ_multiple(:, ismember(curr_succ_multiple(2,:),0)) = [];
    
    % If there are gaps in slices, identify contiguous blocks
    if any(diff(curr_succ_multiple(2,:)) > 1)
        slice_num                     = min(curr_succ_multiple(2,:)) : max(curr_succ_multiple(2,:));
        non_contigous                  = ismember(slice_num, curr_succ_multiple(2,:));
        split_array_consec_int{ii}    = zeros_function_locater(non_contigous); % For Pixel_Index_Per_Slice separations
        split_array_slice{ii}         = cellfun(@(x) slice_num(x), split_array_consec_int{ii}, 'uni', 0); % For Linear_Idx and Final_Centroid separations
    end
    
    % Store fiber matches
    curr_succ_multiple_cell{ii} = curr_succ_multiple;
end
toc

end

%% =============================================================================
% Helper function: zeros_function_locater
% =============================================================================
% Finds contiguous blocks of ones in a logical array
% Input: non_rogue_area_find - logical array indicating presence of fiber in slices
% Output: blocks - cell array of contiguous index ranges
%% =============================================================================
function [blocks] = zeros_function_locater(non_rogue_area_find)

% Multiply by index to preserve positions
rogue_area_find_idx = double(non_rogue_area_find) .* (1:numel(non_rogue_area_find));

% Pad with zeros at ends
wrap = [0, rogue_area_find_idx, 0];

% Diff identifies start and end of contiguous blocks
temp       = diff(wrap ~= 0);
blockStart = find(temp == 1) + 1;
blockEnd   = find(temp == -1);

% Extract contiguous blocks as cell array
blocks = arrayfun(@(bId) wrap(blockStart(bId):blockEnd(bId)), 1:numel(blockStart), 'UniformOutput', false);

end