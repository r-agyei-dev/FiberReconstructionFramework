% INTERSECTOR_CORRECTOR
% -------------------------------------------------------------------------
% PURPOSE:
% Adjusts overlapping ellipses in segmented images by scaling conflicting
% ellipse dimensions and recomputing statistics. Ensures that intersecting
% regions are corrected before final area and pixel calculations.
%
% INPUTS (via varargin):
%   1) STATS_UNCORRECTED  - Original structure of segmented ellipses
%   2) Clean_matrix       - Binary image representing the segmented regions
%   3) C_value            - Scaling factor for correcting intersecting ellipses
%   4) columnsInImage     - Meshgrid of column coordinates (for compatibility)
%   5) rowsInImage        - Meshgrid of row coordinates (for compatibility)
%
% OUTPUT:
%   STATS_CORRECTED_NEW_2 - Corrected statistics structure with updated
%                           ellipse parameters, pixel indices, and areas
% -------------------------------------------------------------------------

function [STATS_CORRECTED_NEW_2] = INTERSECTOR_CORRECTOR(varargin)

STATS_UNCORRECTED  =  varargin{1};
Clean_matrix       =  varargin{2};
C_value            =  varargin{3};
columnsInImage     =  varargin{4}; %#ok<NASGU>
rowsInImage        =  varargin{5}; %#ok<NASGU>

%% ------------------------------------------------------------------------
% Extract Pixel Indices of all ellipses
% -------------------------------------------------------------------------
Pix_idx_corrected = {STATS_UNCORRECTED.Pixel_IDX_List};

%% ------------------------------------------------------------------------
% Create labeled image regions to identify overlaps
% -------------------------------------------------------------------------
Temp_corrected = zeros(size(Clean_matrix));
intersect_checker = cell(1,length(Pix_idx_corrected));

for tt = 1:length(Pix_idx_corrected)
    if isempty(nonzeros(Temp_corrected(Pix_idx_corrected{tt})))
        % No overlap: assign label
        Temp_corrected(Pix_idx_corrected{tt}) = tt;
        intersect_checker{tt} = tt;
    else
        % Overlap detected: record all intersecting labels
        intersect_checker{tt} = vertcat(unique(nonzeros(Temp_corrected(Pix_idx_corrected{tt}))), tt);
    end
end

%% ------------------------------------------------------------------------
% Identify multi-overlap ellipses (more than 1 label)
% and prepare for correction
% -------------------------------------------------------------------------
intersect_checker(cell2mat(cellfun(@(x)length(x)<=1, intersect_checker,'uni',0))) = [];
intersect_checker_idx = vertcat(intersect_checker{:});

%% ------------------------------------------------------------------------
% Apply correction to intersecting ellipses
% Scale their major and minor axes by C_value
% -------------------------------------------------------------------------
STATS_CORRECTED = STATS_UNCORRECTED;

for ii = 1:length(intersect_checker_idx)
    STATS_CORRECTED(intersect_checker_idx(ii)).MajorAxisLength = ...
        C_value * STATS_UNCORRECTED(intersect_checker_idx(ii)).MajorAxisLength;
    STATS_CORRECTED(intersect_checker_idx(ii)).MinorAxisLength = ...
        C_value * STATS_UNCORRECTED(intersect_checker_idx(ii)).MinorAxisLength;
end

%% ------------------------------------------------------------------------
% Recompute final statistics, areas, and pixel indices
% -------------------------------------------------------------------------
[STATS_CORRECTED_NEW_2] = FINAL_STATS_FORMULATOR_SMART_INTERSECT_CORRECTOR( ...
    STATS_UNCORRECTED, STATS_CORRECTED, Clean_matrix, columnsInImage, rowsInImage);

end