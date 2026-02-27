%% =============================================================================
% EPAD_FUNCTION_CODE_Fiber_Stitch_verSPLIT
% =============================================================================
% Description:
% This function identifies overlaps between fibers in consecutive slices,
% separates non-contiguous segments, and prepares arrays for stitching or
% further processing.
%
% INPUTS:
% -------
% curr_vol : 3D array or linear indices representing fibers in the current slice
% succ_vol : 3D array or linear indices representing fibers in the successive slice
%
% OUTPUTS:
% --------
% curr_succ_multiple       : Cell array of [current, successive] fiber matches per fiber
% split_array_consec_int   : Cell array indicating contiguous slice segments
% split_array_slice        : Cell array of slice indices for each contiguous segment
%
% NOTES:
% ------
% - Uses quotient/remainder trick to encode fiber pairs uniquely
% - Handles gaps in fiber slices and splits them into separate contiguous blocks
% - Similar in logic to EPAD_FUNCTION_CODE_Fiber_Stitch_ver
%% =============================================================================

function [curr_succ_multiple, split_array_consec_int, split_array_slice] = EPAD_FUNCTION_CODE_Fiber_Stitch_verSPLIT(curr_vol, succ_vol)

tic

%% ============================================================================
% STEP 1: Encode fiber overlaps using a large exponent
% ============================================================================
% Compute exponent to ensure unique encoding of fiber pairs
exponent = 10^(numel(num2str(max(succ_vol(:)))) + 2);

% Multiply current slice indices to prepare for unique encoding
curr_vol_pre = curr_vol * exponent;

%% ============================================================================
% STEP 2: Combine current and successive slice indices
% ============================================================================
curr_succ_SUM = curr_vol_pre + succ_vol; % Encode pairwise matches

% Clear memory of input volumes to save RAM
disp('clear to relieve memory')
clearvars curr_vol succ_vol

% Find unique non-zero entries representing fiber overlaps
UNIQUE_SUM = unique(nonzeros(curr_succ_SUM(:)));
disp('clear to relieve memory')
clearvars curr_vol_pre curr_succ_SUM

% Remove invalid entries (less than exponent)
UNIQUE_SUM(UNIQUE_SUM < exponent) = [];

%% ============================================================================
% STEP 3: Decode the fiber pairs
% ============================================================================
% Quotient gives the fiber in the current slice
% Remainder gives the corresponding fiber in the successive slice
for ii = 1:numel(UNIQUE_SUM)
    quotient(ii)  = floor(UNIQUE_SUM(ii)/exponent);
    remainder(ii) = mod(UNIQUE_SUM(ii), exponent);
end

% Combine into a two-column array: [current fiber, successive fiber]
curr_succ = [quotient' remainder'];
toc
disp('clear to relieve memory')
clearvars quotient remainder

%% ============================================================================
% STEP 4: Initialize output containers
% ============================================================================
unique_fib_idx = nonzeros(unique(curr_succ(:,1)));      % Unique fibers in current slice
curr_succ_multiple = cell(1, numel(unique_fib_idx));   % Fiber matches per fiber
split_array_consec_int = cell(1, numel(unique_fib_idx)); % Contiguous slice blocks
split_array_slice = cell(1, numel(unique_fib_idx));    % Slice indices for each block

%% ============================================================================
% STEP 5: Extract fiber matches and identify contiguous segments
% ============================================================================
tic
for ii = 1:numel(unique_fib_idx)
    
    % Display progress
    sprintf('%s%d','Extracting curr_succ for fiber ', ii)
    
    % Extract all matches for this fiber
    curr_succ_multiple{ii} = (curr_succ(ismember(curr_succ(:,1), unique_fib_idx(ii)), :))';
    
    % Remove zero entries in successive slice column
    curr_succ_multiple{ii}(:, ismember(curr_succ_multiple{ii}(2,:), 0)) = [];
    
    % Check for non-contiguous slices
    if any(diff(curr_succ_multiple{ii}(2,:)) > 1)
        % Create a full slice range from min to max
        slice_num = min(curr_succ_multiple{ii}(2,:)) : max(curr_succ_multiple{ii}(2,:));
        
        % Identify which slices are non-contiguous
        non_contigous = ismember(slice_num, curr_succ_multiple{ii}(2,:));
        
        % Get contiguous blocks for Pixel_Index_Per_Slice separations
        split_array_consec_int{ii} = zeros_function_locater(non_contigous);
        
        % Get actual slice indices for Linear_Idx and Final_Centroid separations
        split_array_slice{ii} = cellfun(@(x) slice_num(x), split_array_consec_int{ii}, 'uni', 0);
    end
    
end
toc

end

%% =============================================================================
% Helper function: zeros_function_locater
% =============================================================================
% Description:
% Given a logical array, identifies contiguous blocks of ones. Useful for
% splitting fibers that are separated by gaps in slice space.
%
% INPUT:
% non_rogue_area_find : logical array where 1 indicates fiber present
%
% OUTPUT:
% blocks : cell array of contiguous index ranges
%% =============================================================================
function [blocks] = zeros_function_locater(non_rogue_area_find)

% Multiply logical array by indices to preserve positions
rogue_area_find_idx = double(non_rogue_area_find) .* (1:numel(non_rogue_area_find));

% Pad with zeros at both ends to detect edges
wrap = [0, rogue_area_find_idx, 0];

% Use diff to find start (1) and end (-1) of contiguous blocks
temp = diff(wrap ~= 0);
blockStart = find(temp == 1) + 1;
blockEnd   = find(temp == -1);

% Extract contiguous blocks as cell array
blocks = arrayfun(@(bId) wrap(blockStart(bId):blockEnd(bId)), 1:numel(blockStart), 'UniformOutput', false);

end