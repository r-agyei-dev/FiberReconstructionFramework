function [check_mat_center] = SPATIAL_INSTANCE_ASSEMBLER_AUTOMATOR_MULTI_FURTHER(STATS_NEW,iter_slice,size_length)
% SPATIAL_INSTANCE_ASSEMBLER_AUTOMATOR_MULTI_FURTHER
% -------------------------------------------------------------------------
% Converts the pixel coordinates of a fiber (from STATS_NEW) into linear 
% indices within a 3D volume and aggregates them into a single list.
%
% INPUTS:
%   STATS_NEW   : Struct array containing FINAL_VALUES field with pixel coordinates for a fiber
%   iter_slice  : Starting slice index in the 3D volume for this fiber
%   size_length : Size of the 3D volume [rows, cols, slices]
%
% OUTPUT:
%   check_mat_center : Aggregated linear indices for all pixels of this fiber
% -------------------------------------------------------------------------

% Create a temporary copy of the fiber statistics
STATS_NEW_MIDDLE_LINE_temp = STATS_NEW;

hh = iter_slice; % starting slice for this fiber

% Loop over each segment in STATS_NEW and compute linear indices
for jj = 1:size(STATS_NEW_MIDDLE_LINE_temp,2)
    STATS_NEW_MIDDLE_LINE_temp(jj).Pixel_IDX_List = deal( ...
        sub2ind([size_length(1) size_length(2) size_length(3)], ...       % 3D volume size
                STATS_NEW_MIDDLE_LINE_temp(jj).FINAL_VALUES(:,end), ...  % row indices
                STATS_NEW_MIDDLE_LINE_temp(jj).FINAL_VALUES(:,1), ...    % column indices
                (hh-1 + jj) * ones(size(STATS_NEW_MIDDLE_LINE_temp(jj).FINAL_VALUES,1),1))); % slice indices
end

% Aggregate all linear indices from different segments into a single column vector
check_mat_center = vertcat(STATS_NEW_MIDDLE_LINE_temp.Pixel_IDX_List);

% -------------------------------------------------------------------------
% Notes / Debugging Examples:
% STATS_NEW = STATS_NEW_CELL{kk}
% iter_slice = 291
% size_length = [dimX dimY dimZ] of the volume
% The function returns a consolidated list of linear indices for the fiber
% -------------------------------------------------------------------------

end