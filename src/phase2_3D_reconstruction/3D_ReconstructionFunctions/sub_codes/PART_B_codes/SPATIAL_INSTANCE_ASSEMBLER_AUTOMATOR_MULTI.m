function [check_mat_center] = SPATIAL_INSTANCE_ASSEMBLER_AUTOMATOR_MULTI(STATS_NEW_MIDDLE_LINE,fin_lump_VAL_CELL_TEMP_O,iter_slice,size_length)
% SPATIAL_INSTANCE_ASSEMBLER_AUTOMATOR_MULTI
% -------------------------------------------------------------------------
% Converts pixel coordinates of fibers (from STATS_NEW_MIDDLE_LINE) into 
% linear indices within a 3D volume. Handles both single fiber sets and 
% multiple fiber sets stored in a cell array.
%
% INPUTS:
%   STATS_NEW_MIDDLE_LINE          : Struct array or cell array of structs containing FINAL_VALUES (pixel coordinates)
%   fin_lump_VAL_CELL_TEMP_O       : Cell array of offsets for multi-fiber sets
%   iter_slice                     : Starting slice index for this fiber/fiber set
%   size_length                     : Size of the 3D volume [rows, cols, slices]
%
% OUTPUT:
%   check_mat_center               : Linear indices of pixels in the 3D volume
%                                    - Returns a cell array if input is cell
%                                    - Returns a single array if input is a struct array
% -------------------------------------------------------------------------

% Check if input is a cell array of multiple fibers
if iscell(STATS_NEW_MIDDLE_LINE)
    check_mat_center = cell(1,length(STATS_NEW_MIDDLE_LINE));
    
    % Loop through each fiber set in the cell array
    for ii = 1:length(STATS_NEW_MIDDLE_LINE)
        STATS_NEW_MIDDLE_LINE_temp = STATS_NEW_MIDDLE_LINE{ii};
        
        % Compute the slice offset based on the fiber's starting position
        hh = iter_slice + fin_lump_VAL_CELL_TEMP_O{ii}(1) - 1;
        
        % Convert each segment's pixel coordinates to linear indices
        for kk = 1:length(STATS_NEW_MIDDLE_LINE_temp)
            STATS_NEW_MIDDLE_LINE_temp(kk).Pixel_IDX_List = deal( ...
                sub2ind([size_length(1) size_length(2) size_length(3)], ...  % 3D volume dimensions
                        STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES(:,end), ... % row indices
                        STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES(:,1), ...   % column indices
                        (hh-1+kk) * ones(size(STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES,1),1))); % slice indices
        end
        
        % Aggregate all linear indices for this fiber set
        check_mat_center{ii} = vertcat(STATS_NEW_MIDDLE_LINE_temp.Pixel_IDX_List);
    end

else
    % Single fiber set case
    STATS_NEW_MIDDLE_LINE_temp = STATS_NEW_MIDDLE_LINE;
    hh = iter_slice;

    % Convert each segment's pixel coordinates to linear indices
    for kk = 1:size(STATS_NEW_MIDDLE_LINE_temp,2)
        STATS_NEW_MIDDLE_LINE_temp(kk).Pixel_IDX_List = deal( ...
            sub2ind([size_length(1) size_length(2) size_length(3)], ...
                    STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES(:,end), ...
                    STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES(:,1), ...
                    (hh-1+kk) * ones(size(STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES,1),1)));
    end

    % Aggregate all linear indices
    check_mat_center = vertcat(STATS_NEW_MIDDLE_LINE_temp.Pixel_IDX_List);
end