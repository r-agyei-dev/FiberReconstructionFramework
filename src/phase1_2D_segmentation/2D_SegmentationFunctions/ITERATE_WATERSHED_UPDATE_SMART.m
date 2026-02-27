% ITERATE_WATERSHED_UPDATE_SMART
% -------------------------------------------------------------------------
% PURPOSE:
% Iteratively refines segmentation of elliptical regions whose minor axis
% length exceeds a predefined threshold. Performs watershed segmentation
% on oversized regions until they meet the acceptable size criteria.
%
% INPUTS:
%   Image       - Original grayscale or binary image slice
%   level       - Threshold level for watershed segmentation
%   mask_value  - Value used to mask the oversized regions
%   Stats       - Structure array containing ellipse metadata (Centroid,
%                 MajorAxisLength, MinorAxisLength, PixelIdxList, etc.)
%
% OUTPUTS:
%   Stats_keep  - Structure of ellipses that are within size threshold
%   Stats_cell  - Structure of newly segmented regions after watershed
%   Pseudo_Image- Binary image highlighting regions for iterative segmentation
% -------------------------------------------------------------------------

function [Stats_keep,Stats_cell,Pseudo_Image] = ITERATE_WATERSHED_UPDATE_SMART(Image,level, mask_value,Stats)

%% ------------------------------------------------------------------------
% Define minor axis threshold for segmentation refinement
% -------------------------------------------------------------------------
fib_diam_thresh = 10.38;
minor_vec = [Stats.MinorAxisLength];

% Find indices of ellipses exceeding the minor axis threshold
index = find(minor_vec(:) > fib_diam_thresh);

%% ------------------------------------------------------------------------
% Compute area ratios to validate oversized ellipses
% - Compare actual pixel count to ellipse area
% - Keep ellipses that are sufficiently filled
% -------------------------------------------------------------------------
temp = cell2mat(cellfun(@(x)size(x,1), {Stats(index).PixelIdxList}, 'uni', 0)) ./ ...
       (0.25*pi.*[Stats(index).MinorAxisLength].*[Stats(index).MajorAxisLength]);

% Keep ellipses where pixel coverage is at least 90% of ellipse area
temp_idx = index(temp(:) >= 0.90);

% Determine valid and invalid ellipse indices
idx_valid   = sort([find(minor_vec(:) < fib_diam_thresh); temp_idx], 'ascend');
idx_invalid = index(~ismember(index, temp_idx));

%% ------------------------------------------------------------------------
% Separate ellipses that are within acceptable size
% -------------------------------------------------------------------------
Stats_keep = Stats(idx_valid);

%% ------------------------------------------------------------------------
% Generate pseudo-image for iterative watershed refinement of oversized regions
% -------------------------------------------------------------------------
[Pseudo_Image] = PSEUDO_MAT_SMART(Image, Stats, idx_invalid);

%% ------------------------------------------------------------------------
% Apply watershed segmentation on pseudo-image to refine oversized ellipses
% -------------------------------------------------------------------------
[~, Stats_cell] = WATERSHED_FUNC_UPDATE_SMART(Pseudo_Image, level, mask_value);

%% ------------------------------------------------------------------------
% Clean up temporary variables
% -------------------------------------------------------------------------
clearvars minor_vec index temp temp_idx idx_valid idx_invalid
end