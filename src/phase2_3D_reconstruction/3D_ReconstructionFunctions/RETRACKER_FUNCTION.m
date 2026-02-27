%% =========================================================================
% PURPOSE:
% This function searches through a cell array to locate the slice that
% contains a region with a specified size (length_check).
%
% Once the matching slice is found:
%   • slice_region → index/indices where the match occurs
%   • ii_idx       → corresponding tracker/index value from row 2
%
% INPUTS:
%   numel_ellipses → 2×N cell array:
%                     row 1 = region sizes per slice
%                     row 2 = corresponding index mapping
%
%   length_check   → target region size to search for
%
% OUTPUTS:
%   ii_idx        → index value associated with the found slice
%   slice_region  → location(s) within the slice that match length_check
%% =========================================================================

function [ii_idx , slice_region] = RETRACKER_FUNCTION(numel_ellipses,length_check)

% -------------------------------------------------------------------------
% Scan through each slice entry to find the target region size
% -------------------------------------------------------------------------
for ii = 1:length(numel_ellipses)

    % Find positions where the region size equals the target value
    slice_region = find(numel_ellipses{1,ii} == length_check);

    % Stop searching once a match is found
    if ~isempty(slice_region)
        break
    end

end

% -------------------------------------------------------------------------
% Retrieve the corresponding index value from row 2
% -------------------------------------------------------------------------
ii_idx = numel_ellipses{2,ii};