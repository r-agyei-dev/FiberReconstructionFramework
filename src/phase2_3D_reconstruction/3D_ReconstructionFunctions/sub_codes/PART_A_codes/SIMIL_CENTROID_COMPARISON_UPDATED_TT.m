%%%%%%%%%%% THIS FUNCTION COMPUTES SIMILARITY BETWEEN CENTROIDS ACROSS SLICES
% It outputs two key variables:
%   1. sieve_parameter_unique_CELL  -> tracks the indices of matching ellipses across successive slices
%   2. remnant_ellipses             -> indices of unmatched ellipses on each slice

% Brief overview of the steps in this function:
% (1)  Compute distances between centroids of ellipses in consecutive slices using knnsearch.
% (2)  Construct a matrix with each row containing the current ellipse, its distance to candidates, and candidate indices.
% (3)  Identify candidate ellipses with minimal distance.
% (4)  Compare candidate centroids based on minimal X and Y percentage change (optional in further refinement).
% (5)  Resolve conflicts where multiple current ellipses map to the same candidate (competing current ellipses).
% (6)  Set competing or conflicting matches to zero.
% (7)  Store the processed results in sieve_parameter_unique_CELL.
% (8)  Determine remnant (unmatched) ellipses for the next slice.

function [sieve_parameter_unique_CELL, remnant_ellipses] = SIMIL_CENTROID_COMPARISON_UPDATED_TT(varargin)

% Input variables:
% numberOfImageFiles   = number of slices in the 3D image volume
% SLICE_UPDATE         = structure array containing ellipse metadata for each slice
% CENTROID_UPDATE      = cell array containing centroid coordinates of ellipses for each slice

numberOfImageFiles   = varargin{1};
SLICE_UPDATE         = varargin{2};
CENTROID_UPDATE      = varargin{3};

% Initialize output cells
sieve_parameter_unique_CELL = cell(1, size(SLICE_UPDATE, 2));
remnant_ellipses = cell(1, size(SLICE_UPDATE, 2));

% Loop over slices to compute centroid matching (stop two slices before the end)
for tt = 1 : size(SLICE_UPDATE, 2)-2 

    % Clear loop-specific temporary variables
    clearvars dist_matrix dist_matrix_sort logical_idx cord_idx tmp_idx Cen_info_next Cen_info_current diff_array column_percentage_change
    clearvars column_percent_change_cell row_percent_change_cell total_percent_change total_percent_change_cell tmp_idx_vector
    clearvars tmp_idx_vector nt bint multiplet tmp_sort total_percent_change_columns min_val min_val_idx 
    clearvars adel adel_vec sieve_parameter_unique sieve_parameter unique_regions

    current = tt;
    next    = tt + 1;

    % Get number of ellipses on the current and next slices
    length_next    = size(SLICE_UPDATE{next}, 1);       % next slice
    length_current = size(SLICE_UPDATE{current}, 1);    % current slice

    % ===>>> Compute distances between all centroids in current vs next slice
    [Idx, D] = knnsearch(CENTROID_UPDATE{next}, CENTROID_UPDATE{current});

    %%
    % Create sieve_parameter matrix:
    % Column 1: current ellipse number
    % Column 2: distance to candidate ellipse
    % Column 3: index of candidate ellipse
    sieve_parameter = [(1:length_current)' D Idx];

    % Identify competing candidate ellipses (same candidate assigned to multiple current ellipses)
    [nt, ~]          = histc(sieve_parameter(:,3), unique(sieve_parameter(:,3)));
    unique_regions    = unique(sieve_parameter(:,3));
    multiplet         = find(nt > 1);

    % Delete "poorly competing" matches by keeping only the minimal distance
    adel = cell(1, length(multiplet));
    for ii = 1:length(multiplet)
        aa_temp = sieve_parameter(:,3) == unique_regions(multiplet(ii));
        [~, aa_temp] = min(sieve_parameter(aa_temp, 2));  % find minimal distance
        aa_temp_check = find(sieve_parameter(:,3) == unique_regions(multiplet(ii)));
        aa_temp_check(aa_temp) = [];
        adel{ii} = aa_temp_check;
    end

    % Concatenate indices to remove
    adel = adel(~cellfun('isempty', adel));
    adel_vec = vertcat(adel{:});

    % Create the final sieve parameter matrix for the current slice
    sieve_parameter_unique = sieve_parameter;
    sieve_parameter_unique([adel_vec], :) = 0;           % Remove conflicting candidates
    sieve_parameter_unique_CELL{tt} = sieve_parameter_unique;

    % Determine unmatched (remnant) ellipses for the next slice
    remnant_ellipses{tt+1} = find(~ismember(1:length_next, sieve_parameter_unique_CELL{tt}(:, end)));

end

% Handle unmatched ellipses for the first slice (if any)
remnant_ellipses{1} = find(sieve_parameter_unique_CELL{1}(:, 1) == 0);