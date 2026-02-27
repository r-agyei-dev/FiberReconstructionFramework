%% =========================================================================
% PURPOSE:
% This function assembles spatial pixel indices and corresponding height
% values for segmented fiber instances across slices.
%
% It supports both:
%   • Cell input (multiple fiber groups)
%   • Non-cell input (single fiber group)
%
% INPUTS:
%   STATS_NEW_MIDDLE_LINE      → structure/cell of ellipse stats
%   fin_lump_VAL_CELL_TEMP_O   → slice offset information
%   iter_slice                 → current slice index
%   size_length                → size of the 3D volume [rows cols slices]
%
% OUTPUTS:
%   varargout{1} → linear pixel indices of fibers
%   varargout{2} → corresponding height (slice) values
%% =========================================================================

function [varargout] = SPATIAL_INSTANCE_ASSEMBLER_AUTOMATOR_MULTI(varargin)

% -------------------------------------------------------------------------
% Unpack inputs
% -------------------------------------------------------------------------
STATS_NEW_MIDDLE_LINE      =   varargin{1};
fin_lump_VAL_CELL_TEMP_O   =   varargin{2};
iter_slice                 =   varargin{3};
size_length                =   varargin{4};

%% ========================================================================
% CASE 1: Input is a cell array (multiple fiber groups)
%% ========================================================================
if iscell(STATS_NEW_MIDDLE_LINE)

     % Preallocate output containers
     check_mat_center = cell(1,length(STATS_NEW_MIDDLE_LINE));
     height_variable  = cell(1,length(STATS_NEW_MIDDLE_LINE));
     
     % Loop through each fiber group
     for ii = 1:length(STATS_NEW_MIDDLE_LINE)

         STATS_NEW_MIDDLE_LINE_temp = STATS_NEW_MIDDLE_LINE{ii};

         % Compute starting slice height for this group
         hh = iter_slice + fin_lump_VAL_CELL_TEMP_O{ii}(1) - 1;

         % Preallocate height info
         Height_Info_slices = cell(1,length(STATS_NEW_MIDDLE_LINE_temp));

         % -----------------------------------------------------------------
         % Process each ellipse within the group
         % -----------------------------------------------------------------
         for kk = 1:length(STATS_NEW_MIDDLE_LINE_temp)

             % Convert (row,col,slice) to linear indices in 3D volume
             [STATS_NEW_MIDDLE_LINE_temp(kk).Pixel_IDX_List] = deal( ...
                 sub2ind([size_length(1) size_length(2) size_length(3)], ...
                 STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES(:,end), ...
                 STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES(:,1), ...
                 (hh-1+kk) * ones(size(STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES,1),1)));

             % Store slice height for each pixel
             Height_Info_slices{kk} = ...
                 (hh-1+kk) * ones(size(STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES,1),1);

         end

         % Concatenate results for this group
         check_mat_center{ii} = vertcat(STATS_NEW_MIDDLE_LINE_temp.Pixel_IDX_List);
         height_variable{ii}  = vertcat(Height_Info_slices{:});

     end

%% ========================================================================
% CASE 2: Input is NOT a cell array (single fiber group)
%% ========================================================================
else

        STATS_NEW_MIDDLE_LINE_temp = STATS_NEW_MIDDLE_LINE;
        hh = iter_slice;

        % Preallocate height info
        Height_Info_slices = cell(1,length(STATS_NEW_MIDDLE_LINE_temp));

        % -----------------------------------------------------------------
        % Process each ellipse
        % -----------------------------------------------------------------
        for kk = 1:size(STATS_NEW_MIDDLE_LINE_temp,2)

            % Convert to linear indices in the 3D volume
            [STATS_NEW_MIDDLE_LINE_temp(kk).Pixel_IDX_List] = deal( ...
                sub2ind([size_length(1) size_length(2) size_length(3)], ...
                STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES(:,end), ...
                STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES(:,1), ...
                (hh-1+kk) * ones(size(STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES,1),1)));

            % Store slice heights
            Height_Info_slices{kk} = ...
                (hh-1+kk) * ones(size(STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES,1),1);

        end

        % Concatenate outputs
        check_mat_center = vertcat(STATS_NEW_MIDDLE_LINE_temp.Pixel_IDX_List);
        height_variable  = vertcat(Height_Info_slices{:});
end

%% ========================================================================
% OUTPUTS
%% ========================================================================
varargout{1}     = check_mat_center;
varargout{end+1} = height_variable;

end