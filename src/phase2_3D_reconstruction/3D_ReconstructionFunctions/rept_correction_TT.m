%% =========================================================================
% BACKGROUND:
% This routine scans slice-wise region metadata to remove duplicated or
% invalid entries caused by overlapping centroid assignments.
%
% The issue typically arises when multiple regions map to the same pixel
% locations due to centroid misplacement during earlier processing.
%
% The function rebuilds a temporary label map per slice and removes any
% region entries that do not actually appear in the reconstructed map.
%
% OUTPUTS:
%   SLICE_UPDATE    → cleaned structure array of region properties
%   CENTROID_UPDATE → corresponding cleaned centroid array
%
% INPUTS:
%   SLICE_UPDATE    → structure array containing ellipse/region metadata
%   CENTROID_UPDATE → centroid coordinates per slice
%   size_length     → size of the 3D image volume
%% =========================================================================

function [SLICE_UPDATE,CENTROID_UPDATE] = rept_correction_TT(varargin)

% -------------------------------------------------------------------------
% INPUT PARSING
% -------------------------------------------------------------------------
SLICE_UPDATE    = varargin{1};   % Region property structures per slice
CENTROID_UPDATE = varargin{2};   % Corresponding centroid coordinates
size_length     = varargin{3};   % Volume dimensions

% -------------------------------------------------------------------------
% Loop through each slice (last index excluded by original design)
% -------------------------------------------------------------------------
for ii = 1 : size(SLICE_UPDATE,2)-1

    % Progress string (generated but not displayed)
    % sprintf('%s%d','engage rept_idx ',ii)

    % ---------------------------------------------------------------------
    % Create temporary 2D label map for the current slice
    % This map helps identify which region indices actually exist
    % ---------------------------------------------------------------------
    TT = zeros([size_length(1) size_length(2)]);

    % Populate the temporary map using pixel index lists
    for tt = 1:numel(SLICE_UPDATE{ii})

        % Assign region label tt to its pixel locations
        TT(SLICE_UPDATE{ii}(tt).Pixel_IDX_List) = tt;

    end

    % ---------------------------------------------------------------------
    % Determine which region indices never appeared in the map
    %
    % Steps:
    %   1. nonzeros(TT) → active region labels in the map
    %   2. unique(...)  → unique valid region IDs
    %   3. ismember(...) → check which original indices are present
    %   4. ~ismember(...) → mark invalid/duplicate entries
    % ---------------------------------------------------------------------
    del_idx = ~ismember( ...
                    1:numel(SLICE_UPDATE{ii}), ...
                    unique(nonzeros(TT)));

    % ---------------------------------------------------------------------
    % Remove invalid region entries and matching centroid rows
    % ---------------------------------------------------------------------
    if ~isempty(del_idx)
        SLICE_UPDATE{ii}(del_idx)      = [];
        CENTROID_UPDATE{ii}(del_idx,:) = [];
    end

end

end