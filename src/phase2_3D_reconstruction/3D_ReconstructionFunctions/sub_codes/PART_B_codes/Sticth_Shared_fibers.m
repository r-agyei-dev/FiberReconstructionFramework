% STICHTH_SHARED_FIBERS
% This function merges overlapping fibers by averaging their centroid positions
% on shared slices. Metadata of the secondary (probe) fiber is updated to 
% remove the merged slices, and its z-range is recalculated.
%
% Inputs:
%   Centroid_MID_cell_LUMP_ORIG_UP  - Cell array of fiber centroids
%   z_UP                            - Cell array of z-coordinate ranges for each fiber
%   c_f                             - Index of the main (seed) fiber
%   probe_var_over_segTMP           - Index of the fiber to merge with
%   STATS_Pix_LOC_LUMP_UP           - Pixel-level statistics for each fiber

function [varargout] = Sticth_Shared_fibers(varargin)

Centroid_MID_cell_LUMP_ORIG_UP = varargin{1};
z_UP                           = varargin{2};
c_f                            = varargin{3};
probe_var_over_segTMP          = varargin{4};
STATS_Pix_LOC_LUMP_UP          = varargin{5};

% Extract centroid points for seed fiber (p1) and probe fiber (p2)
p1 = Centroid_MID_cell_LUMP_ORIG_UP{c_f};
p2 = Centroid_MID_cell_LUMP_ORIG_UP{probe_var_over_segTMP};

%% Average centroid positions on overlapping slices
% Identify common slices and take element-wise average of centroid coordinates
p1(ismember(p1(:,3), p2(:,3)), :) = 0.5*(p1(ismember(p1(:,3), p2(:,3)), :) + p2(ismember(p2(:,3), p1(:,3)), :));

% Remove merged slices from probe fiber metadata
Centroid_MID_cell_LUMP_ORIG_UP{probe_var_over_segTMP}(ismember(p2(:,3), p1(:,3)), :) = [];
STATS_Pix_LOC_LUMP_UP{probe_var_over_segTMP}(ismember(p2(:,3), p1(:,3))) = [];

% Update the z-range of the probe fiber after removal of merged slices
z_UP{probe_var_over_segTMP} = [min(Centroid_MID_cell_LUMP_ORIG_UP{probe_var_over_segTMP}(:,3)), ...
                                max(Centroid_MID_cell_LUMP_ORIG_UP{probe_var_over_segTMP}(:,3))];

% Return updated centroid of the seed fiber and updated probe fiber metadata
NEW_CENTROID = Centroid_MID_cell_LUMP_ORIG_UP{c_f};

varargout{1} = NEW_CENTROID;
varargout{end+1} = Centroid_MID_cell_LUMP_ORIG_UP{probe_var_over_segTMP};
varargout{end+1} = z_UP{probe_var_over_segTMP};
varargout{end+1} = STATS_Pix_LOC_LUMP_UP{probe_var_over_segTMP};