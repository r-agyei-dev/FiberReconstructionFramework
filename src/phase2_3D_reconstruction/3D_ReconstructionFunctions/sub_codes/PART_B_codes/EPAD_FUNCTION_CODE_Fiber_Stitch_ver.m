function [curr_succ_multiple, split_array_consec_int, split_array_slice] = EPAD_FUNCTION_CODE_Fiber_Stitch_ver(curr_vol, succ_vol)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EPAD_FUNCTION_CODE_Fiber_Stitch_ver
%
% This function identifies overlapping fibers between two consecutive 3D volumes
% (curr_vol and succ_vol) and organizes the overlapping indices for stitching.
%
% Outputs:
% - curr_succ_multiple    : Cell array of overlapping fiber index pairs per fiber
% - split_array_consec_int: Cell array indicating consecutive index blocks for each fiber
% - split_array_slice     : Cell array of slice numbers corresponding to each block
%
% Approach:
% 1. Combine the current and successive volume indices using a large exponent
%    to avoid collisions.
% 2. Extract quotient and remainder to determine the fiber matches.
% 3. Organize overlapping fibers into cell arrays and identify consecutive slices.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tic
% Compute exponent to encode fibers uniquely for summation
exponent = 10^(numel(num2str(max(succ_vol(:)))) + 2);

% Encode and sum the volumes to identify overlaps
curr_vol_pre = curr_vol * exponent;
curr_succ_SUM = curr_vol_pre + succ_vol;

% Clear inputs to save memory
disp('clear to relieve memory')
clearvars curr_vol succ_vol

% Find all unique overlapping sums
UNIQUE_SUM = unique(nonzeros(curr_succ_SUM(:)));
disp('clear to relieve memory')
clearvars curr_vol_pre curr_succ_SUM

% Remove sums that are smaller than exponent (i.e., not overlapping fibers)
UNIQUE_SUM(UNIQUE_SUM < exponent) = [];

% Decode quotient (current volume fiber) and remainder (successor fiber)
for ii = 1:numel(UNIQUE_SUM)
    quotient(ii) = floor(UNIQUE_SUM(ii)/exponent);
    remainder(ii) = mod(UNIQUE_SUM(ii), exponent);
end
curr_succ = [quotient' remainder'];
toc

disp('clear to relieve memory')
clearvars quotient remainder

% Identify all unique fibers in the current volume
unique_fib_idx = nonzeros(unique(curr_succ(:,1)));

% Initialize cell arrays to store results
curr_succ_multiple = cell(1, numel(unique_fib_idx));
split_array_consec_int = cell(1, numel(unique_fib_idx));
split_array_slice = cell(1, numel(unique_fib_idx));

tic
for ii = 1:numel(unique_fib_idx)
    sprintf('%s%d','Extracting curr_succ for fiber ', ii)

    % Extract all matches for the current fiber
    curr_succ_multiple{ii} = curr_succ(ismember(curr_succ(:,1), unique_fib_idx(ii)), :)';
    % Remove zero entries in the successor indices
    curr_succ_multiple{ii}(:, ismember(curr_succ_multiple{ii}(2,:), 0)) = [];

    % Identify non-consecutive slices for separation
    if any(diff(curr_succ_multiple{ii}(2,:)) > 1)
        slice_num = min(curr_succ_multiple{ii}(2,:)) : max(curr_succ_multiple{ii}(2,:));
        non_contigous = ismember(slice_num, curr_succ_multiple{ii}(2,:));
        split_array_consec_int{ii} = zeros_function_locater(non_contigous);                % identify consecutive blocks
        split_array_slice{ii} = cellfun(@(x) slice_num(x), split_array_consec_int{ii}, 'uni', 0); % map blocks to slice numbers
    end
end
toc

end


%% Helper function: locate consecutive nonzero blocks
function [blocks] = zeros_function_locater(non_rogue_area_find)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ZEROS_FUNCTION_LOCATER
%
% Given a logical vector, identifies contiguous nonzero regions and returns them
% as cell array of indices. Used to split fibers across slices.
%
% Input:
% - non_rogue_area_find : logical vector, true where fiber exists
% Output:
% - blocks              : cell array of index vectors corresponding to consecutive regions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Convert logical to actual indices (0 where false)
rogue_area_find_idx = double(non_rogue_area_find) .* (1:numel(non_rogue_area_find));

% Pad with zeros to detect boundaries at ends
wrap = [0, rogue_area_find_idx, 0];

% Identify start and end of each contiguous block
temp = diff(wrap ~= 0);
blockStart = find(temp == 1) + 1;
blockEnd = find(temp == -1);

% Store each block as a cell
blocks = arrayfun(@(bId) wrap(blockStart(bId):blockEnd(bId)), 1:numel(blockStart), 'UniformOutput', false);

end