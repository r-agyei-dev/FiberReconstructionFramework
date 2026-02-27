%%========================================================================
%% FIBER SEPARATING ALGORITHM (DENSE FORMAT, UPDATED)
%%========================================================================
% Version Updates:
% - 2/20/2018: Initial fiber separation function for splitting lumped fibers
% - 3/27/2018: Updated to use overlapping pixels as the disconnect criteria
%
% This function separates "lumped" or overlapping fibers into their distinct parts.
% Handles different scenarios: SPLIT, UNSPLIT, COMBINED, SINGLE STRAND fibers.
%
% Inputs (varargin):
% 1. idx_after_first_thresh   : List of fiber indices after first thresholding
% 2. fib_idx_n_slice_loc_cell : Metadata of fiber indices per slice
% 3. properties               : Structure array containing fiber properties per slice
% 4. size_length              : Reference volume size or slice length
%
% Outputs (varargout):
% 1. fin_lump_VAL_CELL        : Final split/unsplit fiber indices per fiber
% 2. GLOBAL_CENT_TOTAL        : Centroid positions for all fibers (split or unsplit)
% 3. GLOBAL_CENT_UNSPLIT_CELL : Centroids of unsplit fibers
% 4. GLOBAL_CENT_SPLIT_CELL   : Centroids of split fibers
% 5. ROGUE_INDX_CELL_TOTAL    : Indices of rogue areas (problematic fiber regions)
% 6. ROGUE_BEGIN_SLICE        : Slice index where the first rogue area occurs
%%========================================================================

function [varargout] = FIBER_SEPERATING_ALGORITHM_CODE_DENSE_FMT_UPDATED(varargin) 

%-----------------------------------
% Parse input arguments
%-----------------------------------
idx_after_first_thresh     =  varargin{1};
fib_idx_n_slice_loc_cell   =  varargin{2};
properties                 =  varargin{3};
size_length                =  varargin{4};

%-----------------------------------
% Preallocate cell arrays for outputs
%-----------------------------------
GLOBAL_CENT_SPLIT_CELL     = cell(1,length(idx_after_first_thresh));   % Centroids for split fibers
GLOBAL_CENT_UNSPLIT_CELL   = cell(1,length(idx_after_first_thresh));   % Centroids for unsplit fibers
PIX_DIM_CELL               = cell(1,length(idx_after_first_thresh));   % Fiber lengths/dimensions
fin_lump_VAL_CELL          = cell(1,length(idx_after_first_thresh));   % Final indices for split/unsplit fibers
ROGUE_INDX_CELL_S          = cell(1,length(idx_after_first_thresh));   % Rogue slice indices (split)
ROGUE_INDX_CELL_US         = cell(1,length(idx_after_first_thresh));   % Rogue slice indices (unsplit)
ROGUE_BEGIN_SLICE          = zeros(1,length(idx_after_first_thresh));  % Slice index of first rogue occurrence
DISCONNECTER_THRESHOLD     = 0;                                         % Threshold for fiber disconnect (key tolerance)

%-----------------------------------
% Main loop over all fibers
%-----------------------------------
for mm = 1:length(idx_after_first_thresh)

    %-----------------------------------
    % Find fiber metadata for current fiber
    %-----------------------------------
    [tii, slice_region_tmp] = RETRACKER_FUNCTION(fib_idx_n_slice_loc_cell, idx_after_first_thresh(mm));                                            
    iter_slice      = tii;                  % Slice index where root ellipse occurs
    reg_on_iter_slice = slice_region_tmp;   % Metadata location for this fiber

    %-----------------------------------
    % Determine if fiber is "disconnected" using overlapping pixels
    %-----------------------------------
    clearvars cent_info_diff_check
    cent_info_diff_check = FIBER_CONSTITUENT_DISCONNECT_ESTIMATOR(properties, iter_slice, reg_on_iter_slice);  
    cent_info_dist_check_thresh_LD = find(cent_info_diff_check == DISCONNECTER_THRESHOLD);   % Indices of disconnects

    %-----------------------------------
    % Fiber separation: split vs unsplit
    %-----------------------------------
    clearvars fiber_splits GLOBAL_CENT_SPLIT GLOBAL_CENT_UNSPLIT PIX_DIM
    if ~isempty(cent_info_dist_check_thresh_LD)
        % Split fibers: identify distinct strands
        [GLOBAL_CENT_SPLIT, PIX_DIM, ROGUE_INDX_S] = SPLIT_ANALYZER_GLOBAL_CORRECTOR_REMODELED_FOR_ONE( ...
            cent_info_dist_check_thresh_LD, cent_info_diff_check, properties, iter_slice, reg_on_iter_slice, size_length);  
        GLOBAL_CENT_UNSPLIT = [];
        ROGUE_INDX_US = [];
    else
        % Unsplit fibers: single strand or connected fibers
        [GLOBAL_CENT_UNSPLIT, PIX_DIM, ROGUE_INDX_US] = UN_SPLIT_ANALYZER_GLOBAL_CORRECTOR_REMODELED( ...
            properties, iter_slice, reg_on_iter_slice, size_length);
        GLOBAL_CENT_SPLIT = [];
        ROGUE_INDX_S = [];
    end

    %-----------------------------------
    % Store fiber information in output cells
    %-----------------------------------
    GLOBAL_CENT_SPLIT_CELL{mm}   = GLOBAL_CENT_SPLIT;
    GLOBAL_CENT_UNSPLIT_CELL{mm} = GLOBAL_CENT_UNSPLIT;
    PIX_DIM_CELL{mm}             = PIX_DIM;
    ROGUE_INDX_CELL_S{mm}        = ROGUE_INDX_S;
    ROGUE_INDX_CELL_US{mm}       = ROGUE_INDX_US;
    ROGUE_BEGIN_SLICE(mm)        = iter_slice;

    %-----------------------------------
    % Flag if fiber is split
    %-----------------------------------
    if isempty(GLOBAL_CENT_SPLIT_CELL{mm})
        SPLIT_INDICATOR = 0;
    else
        SPLIT_INDICATOR = 1;
    end

    %-----------------------------------
    % Generate final fiber indices after split/unsplit
    %-----------------------------------
    fin_lump_VAL = FINAL_VALUES_ESTIMATOR_CORRECTED_VERSION(PIX_DIM, SPLIT_INDICATOR);
    fin_lump_VAL_CELL{mm} = fin_lump_VAL;

end

%==========================================================================
% Merge split and unsplit centroids
%==========================================================================
GLOBAL_CENT_TOTAL = GLOBAL_CENT_UNSPLIT_CELL;
ROGUE_INDX_CELL_TOTAL = ROGUE_INDX_CELL_US;

GLOBAL_CENT_TOTAL(~cellfun(@isempty, GLOBAL_CENT_SPLIT_CELL)) = GLOBAL_CENT_SPLIT_CELL(~cellfun(@isempty, GLOBAL_CENT_SPLIT_CELL));
ROGUE_INDX_CELL_TOTAL(~cellfun(@isempty, ROGUE_INDX_CELL_S)) = ROGUE_INDX_CELL_S(~cellfun(@isempty, ROGUE_INDX_CELL_S));

%==========================================================================
% Flatten nested cell arrays for centroids and rogue indices
%==========================================================================
for ii = 1:length(GLOBAL_CENT_TOTAL)
    if iscell(GLOBAL_CENT_TOTAL{ii})
        nested_flag = cellfun(@iscell, GLOBAL_CENT_TOTAL{ii});
        if any(nested_flag)
            GLOBAL_CENT_TOTAL{ii} = [GLOBAL_CENT_TOTAL{ii}{:}];
            ROGUE_INDX_CELL_TOTAL{ii} = [ROGUE_INDX_CELL_TOTAL{ii}{:}];
        end
    end
end

%==========================================================================
% Assign output variables
%==========================================================================
varargout{1} = fin_lump_VAL_CELL;
varargout{2} = GLOBAL_CENT_TOTAL;
varargout{3} = GLOBAL_CENT_UNSPLIT_CELL;
varargout{4} = GLOBAL_CENT_SPLIT_CELL;
varargout{5} = ROGUE_INDX_CELL_TOTAL;
varargout{6} = ROGUE_BEGIN_SLICE;

end


%%

%%========================================================================
%% Function 1: RETRACKER_FUNCTION
%%========================================================================
% Purpose: Given the metadata cell array (numel_ellipses) and a fiber index,
%          find the slice index and region index corresponding to that fiber.
%
% Inputs:
% - numel_ellipses : Cell array storing fiber metadata per slice
% - length_check   : Fiber index to search for
%
% Outputs:
% - ii_idx         : Slice index where the fiber occurs
% - slice_region   : Index of fiber in that slice region
%%========================================================================
function [ii_idx, slice_region] = RETRACKER_FUNCTION(numel_ellipses, length_check)
    for ii = 1:length(numel_ellipses)
        slice_region = find(numel_ellipses{1, ii} == length_check); % Check slice for fiber
        if ~isempty(slice_region)
            break; % Stop once fiber is found
        end
    end
    ii_idx = numel_ellipses{2, ii}; % Return slice index
end

%%========================================================================
%% Function 2: FIBER_CONSTITUENT_DISCONNECT_ESTIMATOR
%%========================================================================
% Purpose: Estimate which fibers are "disconnected" based on overlapping pixels.
%          Each consecutive pair of ellipses is checked for shared pixels.
%
% Inputs:
% - properties       : Fiber properties for all slices
% - iter_slice       : Current slice index
% - reg_on_iter_slice: Fiber region index within slice
%
% Output:
% - cent_info_diff_check : Number of overlapping pixels for each consecutive ellipse pair
%%========================================================================
function [cent_info_diff_check] = FIBER_CONSTITUENT_DISCONNECT_ESTIMATOR(properties, iter_slice, reg_on_iter_slice)
    num_fibers = size(properties{iter_slice}{reg_on_iter_slice}, 1);
    cent_info_diff_check = zeros(1, num_fibers - 1);

    for ii = 1:num_fibers-1
        % Count overlapping pixels between consecutive ellipses
        cent_info_diff_check(ii) = numel(nonzeros(ismember( ...
            properties{iter_slice}{reg_on_iter_slice}(ii).Pixel_IDX_List, ...
            properties{iter_slice}{reg_on_iter_slice}(ii+1).Pixel_IDX_List)));
    end
end

%%========================================================================
%% Function 3: SPLIT_ANALYZER_GLOBAL_CORRECTOR_REMODELED_FOR_ONE
%%========================================================================
% Purpose: Analyze split fibers and compute centroids, lengths, and rogue areas.
%
% Inputs:
% - cent_info_dist_check_thresh_LD : Indices of disconnects
% - cent_info_diff_check           : Overlap info for fiber ellipses
% - properties                     : Fiber properties
% - iter_slice                     : Current slice
% - reg_on_iter_slice              : Fiber region index
% - size_length                    : Reference slice length
%
% Outputs:
% - GLOBAL_CENT      : Centroids of split fiber components
% - PIX_DIM          : Lengths of split fiber components
% - ROGUE_INDX_CELL  : Indices of rogue areas
%%========================================================================
function [GLOBAL_CENT, PIX_DIM, ROGUE_INDX_CELL] = SPLIT_ANALYZER_GLOBAL_CORRECTOR_REMODELED_FOR_ONE( ...
        cent_info_dist_check_thresh_LD, cent_info_diff_check, properties, iter_slice, reg_on_iter_slice, size_length)

    n_splits = length(cent_info_dist_check_thresh_LD) + 1;
    fiber_splits = cell(1, n_splits);
    GLOBAL_CENT = cell(1, n_splits);
    PIX_DIM = cell(1, n_splits);
    ROGUE_INDX_CELL = cell(1, n_splits);

    for ii = 1:n_splits
        % Determine indices for this split fiber
        if ii == 1
            fiber_splits{ii} = 1:cent_info_dist_check_thresh_LD(ii);
        elseif ii < n_splits
            fiber_splits{ii} = cent_info_dist_check_thresh_LD(ii-1)+1 : cent_info_dist_check_thresh_LD(ii);
        else
            fiber_splits{ii} = cent_info_dist_check_thresh_LD(ii-1)+1 : length(cent_info_diff_check)+1;
        end

        % Short fiber: just take centroids directly
        if length(fiber_splits{ii}) < 5  % Key tolerance for disconnect
            GLOBAL_CENT{ii} = vertcat(properties{iter_slice}{reg_on_iter_slice}(fiber_splits{ii}).Centroid);
            PIX_DIM{ii} = length(fiber_splits{ii});
            ROGUE_INDX_CELL{ii} = [];
        else
            % Long fiber: invoke contigous band correction for rogue areas
            [PIX_TOT, ~, pix_dim, rogue_indx_find] = CONTIGOUS_BAND_CENTROID_CORRECTOR_REMODELED( ...
                fiber_splits{ii}, properties, iter_slice, reg_on_iter_slice, size_length);
            GLOBAL_CENT{ii} = PIX_TOT;
            PIX_DIM{ii} = pix_dim;
            ROGUE_INDX_CELL{ii} = rogue_indx_find;
        end
    end
end

%%========================================================================
%% Function 4: UN_SPLIT_ANALYZER_GLOBAL_CORRECTOR_REMODELED
%%========================================================================
% Purpose: Analyze a single strand (unsplit) fiber with potential rogue areas.
%
% Inputs/Outputs are similar to SPLIT_ANALYZER but for single strand
%%========================================================================
function [GLOBAL_CENT_UNSPLIT, PIX_DIM, rogue_indx_find] = UN_SPLIT_ANALYZER_GLOBAL_CORRECTOR_REMODELED(properties, iter_slice, reg_on_iter_slice, size_length)
    fiber_splits_vec = 1:length([properties{iter_slice}{reg_on_iter_slice}.MinorAxisLength]);
    [GLOBAL_CENT_UNSPLIT, ~, PIX_DIM, rogue_indx_find] = CONTIGOUS_BAND_CENTROID_CORRECTOR_REMODELED( ...
        fiber_splits_vec, properties, iter_slice, reg_on_iter_slice, size_length);
end

%%========================================================================
%% Function 5: CONTIGOUS_BAND_CENTROID_CORRECTOR_REMODELED
%%========================================================================
% Purpose: Align fiber centroids, interpolate rogue areas, and compute vestiges
%
% Inputs:
% - fiber_splits_vec : indices of fiber ellipses
% - properties       : fiber properties
% - iter_slice       : current slice
% - reg_on_iter_slice: fiber region index
% - size_length      : reference slice length
%
% Outputs:
% - PIX_TOT          : Final centroid coordinates
% - pixels           : Cell of centroid sub-strands
% - pix_dim          : Lengths of individual strands
% - rogue_indx_find  : Indices of rogue areas
%%========================================================================
function [PIX_TOT, pixels, pix_dim, rogue_indx_find] = CONTIGOUS_BAND_CENTROID_CORRECTOR_REMODELED( ...
        fiber_splits_vec, properties, iter_slice, reg_on_iter_slice, size_length)

    % Extract centroid and major axis lengths
    total_centroid_temp = vertcat(properties{iter_slice}{reg_on_iter_slice}(fiber_splits_vec).Centroid);
    total_area_temp = [properties{iter_slice}{reg_on_iter_slice}(fiber_splits_vec).MajorAxisLength];

    % Detect rogue areas based on area distribution
    threshold = 0.09; strict_tmp_val = 0.25; strict_var = 0.5; ease_var = 2.5;
    min_to_mode = mode(round(total_area_temp)) - min(total_area_temp);
    mode_to_max = max(total_area_temp) - mode(round(total_area_temp));
    numerator = min(mode_to_max, min_to_mode);
    denomenator = max(mode_to_max, min_to_mode);

    if (numerator/denomenator) > threshold
        [non_rogue_area_find] = Rogue_region_FUNC(total_area_temp);
    else
        [non_rogue_area_find] = Rogue_region_FUNC_2(total_area_temp, strict_tmp_val, strict_var, ease_var);
    end
    rogue_indx_find = find(~non_rogue_area_find);

    % Interpolate centroids if more than 10 points
    if size(total_centroid_temp,1) > 10
        sample_points_cr = [find(non_rogue_area_find)]';
        sample_values_cr = total_centroid_temp(sample_points_cr, :);
        query_vec_points = 1:size(total_centroid_temp,1);
        result(:,1) = interp1(sample_points_cr, sample_values_cr(:,1), query_vec_points', 'linear','extrap');
        result(:,2) = interp1(sample_points_cr, sample_values_cr(:,2), query_vec_points', 'linear','extrap');
        Final_Result = rogue_extrapolation(non_rogue_area_find, result, total_centroid_temp, size_length);
    else
        Final_Result = total_centroid_temp;
    end

    PIX_TOT = Final_Result;
    pixels = Final_Result;
    pix_dim = size(Final_Result,1);
end

%%========================================================================
%% Function 6: FINAL_VALUES_ESTIMATOR_CORRECTED_VERSION
%%========================================================================
% Purpose: Generate unique indices for fiber strands, accounting for:
% - Split fibers with/without conjoined volumes
% - Unsplit fibers with/without conjoined volumes
%%========================================================================
function [IDX_CONST_ELL] = FINAL_VALUES_ESTIMATOR_CORRECTED_VERSION(PIX_DIM, SPLIT_INDICATOR)
    if ~iscell(PIX_DIM)
        IDX_CONST_ELL = 1:PIX_DIM; % Unsplit single fiber
    elseif iscell(PIX_DIM) && SPLIT_INDICATOR == 1
        temp = cellfun(@iscell, PIX_DIM);
        if ~any(temp)
            % Split fibers without conjoined volumes
            IDX_CONST_ELL = cell(1,length(PIX_DIM));
            for ii = 1:length(PIX_DIM)
                if ii == 1
                    IDX_CONST_ELL{ii} = 1:PIX_DIM{ii};
                else
                    IDX_CONST_ELL{ii} = IDX_CONST_ELL{ii-1}(end)+1 : IDX_CONST_ELL{ii-1}(end) + sum(PIX_DIM{ii});
                end
            end
        else
            % Split fibers with conjoined volumes
            IDX_CONST_ELL = FINAL_VALUES_ESTIMATOR_SUBFILE_SPLIT(PIX_DIM);
        end
    else
        % Unsplit fibers with conjoined volumes
        IDX_CONST_ELL = FINAL_VALUES_ESTIMATOR_SUBFILE_UNSPLIT(PIX_DIM);
    end
end

%%========================================================================
%% Function 7: FINAL_VALUES_ESTIMATOR_SUBFILE_SPLIT
%%========================================================================
% Purpose: Assign unique indices for split fibers, including branched strands
%%========================================================================
function [fin_lump] = FINAL_VALUES_ESTIMATOR_SUBFILE_SPLIT(PIX_DIM)
    if iscell(PIX_DIM)
        PIX_DIM_LUMP_cell = [PIX_DIM{:}];
        TOT_SPECIFIC_LENGTH = length(PIX_DIM_LUMP_cell);
        PIX_DIM_TMP = PIX_DIM;

        extract_ff = cell(1,length(PIX_DIM_TMP));
        extract_bef = cell(1,length(PIX_DIM_TMP));

        for ii = 1:length(PIX_DIM_TMP)
            if iscell(PIX_DIM_TMP{ii})
                for tt = 2:length(PIX_DIM_TMP{ii})
                    extract_ff{ii}(tt-1) = PIX_DIM_TMP{ii}{tt}(1);
                    PIX_DIM_TMP{ii}{tt}(1) = [];
                end
                extract_bef{ii} = [0 -extract_ff{ii}];
            else
                extract_bef{ii} = [];
            end
        end

        if length(find(cellfun(@isempty, extract_bef))) == TOT_SPECIFIC_LENGTH
            fin_lump = zeros(1,TOT_SPECIFIC_LENGTH);
            for ii = 1:TOT_SPECIFIC_LENGTH
                if ii == 1
                    fin_lump{ii} = 1:PIX_DIM_LUMP_cell(ii);
                else
                    fin_lump{ii} = fin_lump{ii-1}(end)+1 : fin_lump{ii-1}(end) + PIX_DIM_LUMP_cell(ii);
                end
            end
        else
            zero_idx = find(cellfun(@isempty, extract_bef));
            for yy = 1:length(zero_idx)
                extract_bef{zero_idx(yy)} = 0;
            end
            extract_bef = [extract_bef{:}];
            fin_lump = cell(1,TOT_SPECIFIC_LENGTH);
            for ii = 1:TOT_SPECIFIC_LENGTH
                if ii ==1
                    fin_lump{ii} = 1:sum(PIX_DIM_LUMP_cell{ii});
                else
                    fin_lump{ii} = fin_lump{ii-1}(end)+extract_bef(ii)+1 : sum([PIX_DIM_LUMP_cell{1:ii}]) + sum(extract_bef(2:ii));
                end
            end
        end
    else
        fin_lump = 1:PIX_DIM;
    end
end

%%========================================================================
%% Function 8: FINAL_VALUES_ESTIMATOR_SUBFILE_UNSPLIT
%%========================================================================
% Purpose: Assign unique indices for unsplit fibers with potential conjoined volumes
%%========================================================================
function [fin_lump] = FINAL_VALUES_ESTIMATOR_SUBFILE_UNSPLIT(PIX_DIM)
    TOT_SPECIFIC_LENGTH = length(PIX_DIM);
    tracker = zeros(1,TOT_SPECIFIC_LENGTH-1);
    for ii = 1:TOT_SPECIFIC_LENGTH-1
        tracker(ii) = PIX_DIM{ii}(end);
    end
    fin_lump = cell(1,TOT_SPECIFIC_LENGTH);
    for ii = 1:TOT_SPECIFIC_LENGTH
        if ii == 1
            fin_lump{ii} = 1:sum(PIX_DIM{ii});
        else
            fin_lump{ii} = fin_lump{ii-1}(end)-tracker(ii-1)+1 : fin_lump{ii-1}(end)-tracker(ii-1) + sum(PIX_DIM{ii});
        end
    end
end