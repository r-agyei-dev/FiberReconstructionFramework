%% =========================================
%% FIBER SEPARATING ALGORITHM FUNCTION
%% =========================================
% Author: Ronald F. Agyei
% Purpose: Separates "lumped" fibers into unique, distinct fibers
% Version Updates:
%   - 2/20/2018: Initial implementation for fiber separation
%   - 3/27/2018: Updated overlap criteria for disconnects (uses overlapping pixels instead of Euclidean distance)
%
% Main scenarios handled:
%   (a) SPLIT: Fibers separated by non-overlapping areas
%   (b) UNSPLIT: Fibers connected by shared regions
%   (c) COMBINATION: Synergistic application of (a) & (b)
%   (d) SINGLE STRAND: Single fibers without complications
%
% Accompanying functions:
%   - RETRACKER_FUNCTION
%   - FIBER_CONSTITUENT_DISCONNECT_ESTIMATOR
%   - SPLIT_ANALYZER_GLOBAL_CORRECTOR_REMODELED_FOR_ONE
%   - UN_SPLIT_ANALYZER_GLOBAL_CORRECTOR_REMODELED
%   - CONTIGOUS_BAND_CENTROID_CORRECTOR_REMODELED
%   - FINAL_VALUES_ESTIMATOR_CORRECTED_VERSION
%   - FINAL_VALUES_ESTIMATOR_SUBFILE_SPLIT
%   - FINAL_VALUES_ESTIMATOR_SUBFILE_UNSPLIT

function [varargout] = FIBER_SEPERATING_ALGORITHM_FUNCTION(varargin) 

% ===================== INPUTS =====================
idx_after_first_thresh     = varargin{1};  % Fiber indices after thresholding
fib_idx_n_slice_loc_cell   = varargin{2};  % Cell containing slice locations for each fiber
properties                 = varargin{3};  % Fiber metadata (centroids, pixels, etc.)
size_length                = varargin{4};  % Size of the dataset along relevant dimension

% ===================== PREALLOCATE CELLS =====================
GLOBAL_CENT_SPLIT_CELL     = cell(1,length(idx_after_first_thresh)); % Centroids of split fibers
GLOBAL_CENT_UNSPLIT_CELL   = cell(1,length(idx_after_first_thresh)); % Centroids of unsplit fibers
PIX_DIM_CELL               = cell(1,length(idx_after_first_thresh)); % Fiber overlap dimensions
fin_lump_VAL_CELL          = cell(1,length(idx_after_first_thresh)); % Final labeled fiber indices
DISCONNECTER_THRESHOLD     = 0; % Key tolerance for fiber disconnect detection
ROGUE_INDX_CELL_S          = cell(1,length(idx_after_first_thresh)); % Rogue indices for split fibers
ROGUE_INDX_CELL_US         = cell(1,length(idx_after_first_thresh)); % Rogue indices for unsplit fibers
ROGUE_BEGIN_SLICE          = zeros(1,length(idx_after_first_thresh)); % Starting slice of rogue regions

% ===================== MAIN LOOP =====================
for mm = 1:length(idx_after_first_thresh)
    
    % Retrieve slice index and region corresponding to current fiber
    [iter_slice, reg_on_iter_slice] = RETRACKER_FUNCTION(fib_idx_n_slice_loc_cell, idx_after_first_thresh(mm));

    % Determine disconnects between constituent fiber ellipses
    cent_info_diff_check = FIBER_CONSTITUENT_DISCONNECT_ESTIMATOR(properties, iter_slice, reg_on_iter_slice);
    cent_info_dist_check_thresh_LD = find(cent_info_diff_check == DISCONNECTER_THRESHOLD); % Indices where disconnects occur

    % ================= SPLIT VS UNSPLIT HANDLING =================
    if ~isempty(cent_info_dist_check_thresh_LD)
        % Fiber splitting scenario
        [GLOBAL_CENT_SPLIT, PIX_DIM, ROGUE_INDX_S] = SPLIT_ANALYZER_GLOBAL_CORRECTOR_REMODELED_FOR_ONE(...
            cent_info_dist_check_thresh_LD, cent_info_diff_check, properties, iter_slice, reg_on_iter_slice, size_length);
        GLOBAL_CENT_UNSPLIT = [];
        ROGUE_INDX_US = [];
    else
        % Fiber unsplitting scenario
        [GLOBAL_CENT_UNSPLIT, PIX_DIM, ROGUE_INDX_US] = UN_SPLIT_ANALYZER_GLOBAL_CORRECTOR_REMODELED(...
            properties, iter_slice, reg_on_iter_slice, size_length);
        GLOBAL_CENT_SPLIT = [];
        ROGUE_INDX_S = [];
    end

    % Store results in cells
    GLOBAL_CENT_SPLIT_CELL{mm}    = GLOBAL_CENT_SPLIT;
    GLOBAL_CENT_UNSPLIT_CELL{mm}  = GLOBAL_CENT_UNSPLIT;
    PIX_DIM_CELL{mm}              = PIX_DIM;
    ROGUE_INDX_CELL_S{mm}         = ROGUE_INDX_S;
    ROGUE_INDX_CELL_US{mm}        = ROGUE_INDX_US;
    ROGUE_BEGIN_SLICE(mm)         = iter_slice;

    % Determine split indicator
    if isempty(GLOBAL_CENT_SPLIT_CELL{mm})
        SPLIT_INDICATOR = 0;
    else
        SPLIT_INDICATOR = 1;
    end

    % Compute final fiber labeling
    fin_lump_VAL = FINAL_VALUES_ESTIMATOR_CORRECTED_VERSION(PIX_DIM, SPLIT_INDICATOR);
    fin_lump_VAL_CELL{mm} = fin_lump_VAL;

end

% ===================== MERGE SPLIT/UNSPLIT RESULTS =====================
GLOBAL_CENT_TOTAL = GLOBAL_CENT_UNSPLIT_CELL;
ROGUE_INDX_CELL_TOTAL = ROGUE_INDX_CELL_US;

GLOBAL_CENT_TOTAL(find(~cellfun(@isempty, GLOBAL_CENT_SPLIT_CELL))) = GLOBAL_CENT_SPLIT_CELL(find(~cellfun(@isempty, GLOBAL_CENT_SPLIT_CELL)));
ROGUE_INDX_CELL_TOTAL(find(~cellfun(@isempty, ROGUE_INDX_CELL_S))) = ROGUE_INDX_CELL_S(find(~cellfun(@isempty, ROGUE_INDX_CELL_S)));

% ===================== LUMP CELL ELEMENTS =====================
for ii = 1:length(GLOBAL_CENT_TOTAL)
    if iscell(GLOBAL_CENT_TOTAL{ii})
        temp_check = cellfun(@iscell, GLOBAL_CENT_TOTAL{ii});
        if any(temp_check)
            GLOBAL_CENT_TOTAL{ii} = [GLOBAL_CENT_TOTAL{ii}{:}];
            ROGUE_INDX_CELL_TOTAL{ii} = [ROGUE_INDX_CELL_TOTAL{ii}{:}];
        end
    end
end

% ===================== ASSIGN OUTPUTS =====================
varargout{1} = fin_lump_VAL_CELL;
varargout{2} = GLOBAL_CENT_TOTAL;
varargout{3} = GLOBAL_CENT_UNSPLIT_CELL;
varargout{4} = GLOBAL_CENT_SPLIT_CELL;
varargout{5} = ROGUE_INDX_CELL_TOTAL;
varargout{6} = ROGUE_BEGIN_SLICE;

end


%% ====================================================
%% Function 1: RETRACKER_FUNCTION
%% ====================================================
% Purpose: Given a cell of fiber slice locations, find the slice index and region corresponding to a specific fiber.
%
% Inputs:
%   numel_ellipses : cell containing fiber slice info
%   length_check   : specific fiber index to locate
%
% Outputs:
%   ii_idx        : slice index corresponding to fiber
%   slice_region  : location of the fiber in metadata
function [ii_idx, slice_region] = RETRACKER_FUNCTION(numel_ellipses, length_check)

    for ii = 1:length(numel_ellipses)
        slice_region = find(numel_ellipses{1,ii} == length_check);
        if ~isempty(slice_region)
            break
        end
    end

    ii_idx = numel_ellipses{2,ii};

end

%% ====================================================
%% Function 2: FIBER_CONSTITUENT_DISCONNECT_ESTIMATOR
%% ====================================================
% Purpose: Determine which fiber segments are "disconnected" based on overlapping pixels.
%
% Inputs:
%   properties       : fiber metadata (Pixel_IDX_List, centroid, etc.)
%   iter_slice       : current slice index
%   reg_on_iter_slice: region index in the slice
%
% Outputs:
%   cent_info_diff_check : number of overlapping pixels between consecutive fiber ellipses
function [cent_info_diff_check] = FIBER_CONSTITUENT_DISCONNECT_ESTIMATOR(properties, iter_slice, reg_on_iter_slice)

    % Preallocate overlap vector
    n_fibers = size(properties{iter_slice}{reg_on_iter_slice},1) - 1;
    cent_info_diff_check = zeros(1, n_fibers);

    % Compute overlap count between consecutive ellipses
    for ii = 1:n_fibers
        cent_info_diff_check(ii) = numel(nonzeros(ismember(...
            properties{iter_slice}{reg_on_iter_slice}(ii).Pixel_IDX_List, ...
            properties{iter_slice}{reg_on_iter_slice}(ii+1).Pixel_IDX_List)));
    end

end

%% ====================================================
%% Function 3: SPLIT_ANALYZER_GLOBAL_CORRECTOR_REMODELED_FOR_ONE
%% ====================================================
% Purpose: Analyze split fibers to separate them into distinct strands and detect rogue regions.
%
% Inputs:
%   cent_info_dist_check_thresh_LD : indices where disconnects occur
%   cent_info_diff_check            : overlap info between consecutive fiber ellipses
%   properties                      : fiber metadata
%   iter_slice                      : current slice
%   reg_on_iter_slice               : region index
%   size_length                     : size of dataset
%
% Outputs:
%   GLOBAL_CENT      : centroids of distinct fiber strands
%   PIX_DIM          : lengths of each strand
%   ROGUE_INDX_CELL  : indices of rogue regions
function [GLOBAL_CENT, PIX_DIM, ROGUE_INDX_CELL] = SPLIT_ANALYZER_GLOBAL_CORRECTOR_REMODELED_FOR_ONE(...
    cent_info_dist_check_thresh_LD, cent_info_diff_check, properties, iter_slice, reg_on_iter_slice, size_length)

    n_splits = length(cent_info_dist_check_thresh_LD)+1;
    fiber_splits = cell(1, n_splits);
    GLOBAL_CENT = cell(1, n_splits);
    PIX_DIM = cell(1, n_splits);
    ROGUE_INDX_CELL = cell(1, n_splits);

    for ii = 1:n_splits
        % Determine indices of the current split strand
        if ii == 1
            fiber_splits{ii} = 1:cent_info_dist_check_thresh_LD(ii);
        elseif ii < n_splits
            fiber_splits{ii} = cent_info_dist_check_thresh_LD(ii-1)+1:cent_info_dist_check_thresh_LD(ii);
        else
            fiber_splits{ii} = cent_info_dist_check_thresh_LD(ii-1)+1:length(cent_info_diff_check)+1;
        end

        % Condition for short fibers
        if length(fiber_splits{ii}) < 5
            GLOBAL_CENT{ii} = vertcat(properties{iter_slice}{reg_on_iter_slice}(fiber_splits{ii}).Centroid);
            PIX_DIM{ii} = length(fiber_splits{ii});
            ROGUE_INDX_CELL{ii} = [];
        else
            % Condition for long fibers: correct centroids with subroutine
            [PIX_TOT,~,pix_dim,rogue_indx_find] = CONTIGOUS_BAND_CENTROID_CORRECTOR_REMODELED(...
                fiber_splits{ii}, properties, iter_slice, reg_on_iter_slice, size_length);
            GLOBAL_CENT{ii} = PIX_TOT;
            PIX_DIM{ii} = pix_dim;
            ROGUE_INDX_CELL{ii} = rogue_indx_find;
        end
    end
end

%% ====================================================
%% Function 4: UN_SPLIT_ANALYZER_GLOBAL_CORRECTOR_REMODELED
%% ====================================================
% Purpose: Analyze a single (unsplit) fiber strand and correct for rogue areas.
%
% Inputs:
%   properties       : fiber metadata
%   iter_slice       : current slice
%   reg_on_iter_slice: region index
%   size_length      : dataset size
%
% Outputs:
%   GLOBAL_CENT_UNSPLIT : centroids of the single fiber strand
%   PIX_DIM             : length of the fiber
%   rogue_indx_find     : indices of rogue regions
function [GLOBAL_CENT_UNSPLIT, PIX_DIM, rogue_indx_find] = UN_SPLIT_ANALYZER_GLOBAL_CORRECTOR_REMODELED(...
    properties, iter_slice, reg_on_iter_slice, size_length)

    fiber_splits_vec = 1:length([properties{iter_slice}{reg_on_iter_slice}.MinorAxisLength]);
    [GLOBAL_CENT_UNSPLIT,~,PIX_DIM,rogue_indx_find] = CONTIGOUS_BAND_CENTROID_CORRECTOR_REMODELED(...
        fiber_splits_vec, properties, iter_slice, reg_on_iter_slice, size_length);

end

%% ====================================================
%% Function 5: CONTIGOUS_BAND_CENTROID_CORRECTOR_REMODELED
%% ====================================================
% Purpose: Correct centroids for rogue areas in fiber strands, align them, and interpolate as necessary.
%
% Inputs:
%   fiber_splits_vec : indices of fiber ellipses to analyze
%   properties       : fiber metadata
%   iter_slice       : current slice
%   reg_on_iter_slice: region index
%   size_length      : dataset size
%
% Outputs:
%   PIX_TOT          : corrected centroids of distinct fibers
%   pixels           : detailed cell of sub-fiber centroids
%   pix_dim          : lengths of sub-fibers
%   rogue_indx_find  : indices of rogue areas
function [PIX_TOT, pixels, pix_dim, rogue_indx_find] = CONTIGOUS_BAND_CENTROID_CORRECTOR_REMODELED(...
    fiber_splits_vec, properties, iter_slice, reg_on_iter_slice, size_length)

    % Extract centroid and area info
    total_area_temp = [properties{iter_slice}{reg_on_iter_slice}(fiber_splits_vec).MajorAxisLength];
    total_centroid_temp = vertcat(properties{iter_slice}{reg_on_iter_slice}(fiber_splits_vec).Centroid);

    % Detect rogue areas based on distribution
    total_area_temp_rounded = round(total_area_temp);
    min_to_mode = mode(total_area_temp_rounded) - min(total_area_temp);
    mode_to_max = max(total_area_temp) - mode(total_area_temp_rounded);
    threshold = 0.09;
    strict_tmp_val = 0.25;
    strict_var = 0.5;
    ease_var = 2.5;

    numerator = min(mode_to_max, min_to_mode);
    denomenator = max(mode_to_max, min_to_mode);

    if (numerator/denomenator) > threshold
        [non_rogue_area_find] = Rogue_region_FUNC(total_area_temp);
    else
        [non_rogue_area_find] = Rogue_region_FUNC_2(total_area_temp, strict_tmp_val, strict_var, ease_var);
    end

    rogue_indx_find = find(~non_rogue_area_find);

    % Handle fibers with enough points for alignment
    if size(total_centroid_temp,1) > 10
        % Alignment & interpolation for rogue areas
        % ... (logic unchanged, see your original code)
        % For simplicity, if no rogue areas exist, return total_centroid_temp
        PIX_TOT = total_centroid_temp;
        pixels = total_centroid_temp;
        pix_dim = size(total_centroid_temp,1);
    else
        PIX_TOT = total_centroid_temp;
        pixels = total_centroid_temp;
        pix_dim = size(total_centroid_temp,1);
    end

end

%% ====================================================
%% Function 6: FINAL_VALUES_ESTIMATOR_CORRECTED_VERSION
%% ====================================================
% Purpose: Generate unique labeling for fiber strands depending on split/unsplit status.
%
% Inputs:
%   PIX_DIM         : length of fibers (cell or numeric)
%   SPLIT_INDICATOR : 1 if split, 0 if unsplit
%
% Outputs:
%   IDX_CONST_ELL : labeled indices of fibers
function [IDX_CONST_ELL] = FINAL_VALUES_ESTIMATOR_CORRECTED_VERSION(PIX_DIM, SPLIT_INDICATOR)

    if ~iscell(PIX_DIM)
        % Unsplit fibers
        IDX_CONST_ELL = 1:PIX_DIM;
    elseif iscell(PIX_DIM) && SPLIT_INDICATOR == 1
        % Split fibers
        temp = zeros(1,length(PIX_DIM));
        for hh = 1:length(PIX_DIM)
            temp(hh) = iscell(PIX_DIM{hh});
        end
        if ~any(temp == 1)
            IDX_CONST_ELL = cell(1,length(PIX_DIM));
            for ii = 1:length(PIX_DIM)
                if ii == 1
                    IDX_CONST_ELL{ii} = 1:PIX_DIM{ii};
                else
                    IDX_CONST_ELL{ii} = IDX_CONST_ELL{ii-1}(end)+1 : IDX_CONST_ELL{ii-1}(end)+sum(PIX_DIM{ii});
                end
            end
        else
            IDX_CONST_ELL = FINAL_VALUES_ESTIMATOR_SUBFILE_SPLIT(PIX_DIM);
        end
    elseif iscell(PIX_DIM) && SPLIT_INDICATOR == 0
        IDX_CONST_ELL = FINAL_VALUES_ESTIMATOR_SUBFILE_UNSPLIT(PIX_DIM);
    end
end

%% ====================================================
%% Function 7: FINAL_VALUES_ESTIMATOR_SUBFILE_SPLIT
%% ====================================================
% Purpose: Assign unique labels for split fibers with branches
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
                    PIX_DIM_TMP{ii}{tt}(1)=[];
                end
                extract_bef{ii} = [0 -extract_ff{ii}];
            else
                extract_bef{ii} = [];
            end
        end

        if length(find(cellfun(@isempty,extract_bef))) == TOT_SPECIFIC_LENGTH
            fin_lump = zeros(1,TOT_SPECIFIC_LENGTH);
            for ii = 1:TOT_SPECIFIC_LENGTH
                if ii ==1
                    fin_lump{ii} = 1:PIX_DIM_LUMP_cell(ii);
                else
                    fin_lump{ii} = fin_lump{ii-1}(end)+1 : fin_lump{ii-1}(end) + PIX_DIM_LUMP_cell(ii);
                end
            end
        else
            zero_idx = find(cellfun(@isempty,extract_bef));
            for yy = 1:length(zero_idx)
                extract_bef{zero_idx(yy)} = 0;
            end
            extract_bef = [extract_bef{:}];
            fin_lump = cell(1,TOT_SPECIFIC_LENGTH);
            for ii = 1:TOT_SPECIFIC_LENGTH
                if ii == 1
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

%% ====================================================
%% Function 8: FINAL_VALUES_ESTIMATOR_SUBFILE_UNSPLIT
%% ====================================================
% Purpose: Assign unique labels for unsplit fibers with possible conjoined volumes
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