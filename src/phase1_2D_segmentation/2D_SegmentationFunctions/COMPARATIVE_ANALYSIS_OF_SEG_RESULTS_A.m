% This function compares segmentations from sharp and unsharp images and
% selects the most appropriate regions to retain.
%
% The routine classifies regions into four categories:
% (a) Areas captured only by the sharp image
% (b) Areas captured only by the unsharp image
% (c) Competing regions where a preferred candidate is selected
% (d) Ambiguous or unmatched regions handled separately
%
% Design philosophy:
% When sharp and unsharp regions overlap, preference is typically given to
% the more conservative (usually smaller) bounding region because such
% segmentations tend to be less sensitive to blur at ellipse boundaries.

function [varargout] = COMPARATIVE_ANALYSIS_OF_SEG_RESULTS_A(varargin)

% ---------------------- INPUT DEFINITIONS ----------------------
% Stats_Current     : Metadata structure for the current adjusted slice
% Stats_Consecutive : Metadata structure for the consecutive adjusted slice
% Clean_matrix      : Original grayscale image
% Comp_Bin_Im       : Label image storing optimized ellipse segmentation
%                     (updated throughout processing)
% orig_comp         : Flag indicating if this is the first comparison
%                     (1 = initial population, 0 = subsequent updates)
% Base_AREA_compare : Reference areas used to detect over/under segmentation
% columnsInImage    : Meshgrid column indices of the 2D image
% rowsInImage       : Meshgrid row indices of the 2D image

Stats_Current         =    varargin{1};
Stats_Consecutive     =    varargin{2};
Clean_matrix          =    varargin{3};
Comp_Bin_Im           =    varargin{4};
orig_comp             =    varargin{5};
Base_AREA_compare     =    varargin{6};
columnsInImage        =    varargin{7};
rowsInImage           =    varargin{8};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1: Build labeled images for current and consecutive segmentations
% Purpose:
%   Convert pixel index lists into labeled images so overlap operations
%   can be performed efficiently.

base_area_length     =    numel(Base_AREA_compare);

% Extract pixel index lists from stats structures
Pix_idx_curr         =    {Stats_Current.Pixel_IDX_List};
Pix_idx_consec       =    {Stats_Consecutive.Pixel_IDX_List};

% Initialize temporary label images
Temp_sharp_curr            =    zeros(size(Clean_matrix));
Temp_sharp_consec          =    zeros(size(Clean_matrix));

% Populate label image for current slice
for tt = 1:length(Pix_idx_curr)
    Temp_sharp_curr(Pix_idx_curr{tt}) = tt;
end

% Populate label image for consecutive slice
for tt = 1:length(Pix_idx_consec)
    Temp_sharp_consec(Pix_idx_consec{tt}) = tt;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2: Extract region areas and compute overlaps

% Area vectors for each segmentation set
areas_sharp_curr   = [Stats_Current.Final_Area];
areas_sharp_consec = [Stats_Consecutive.Final_Area];

% Determine overlapping and unmatched regions between slices
[Curr_Cell, Curr_uncaptured_cell, Consec_uncaptured_cell, unsure_probe_vals] = ...
    EPAD_FUNCTION_CODE_MASTER_VERSION_2(Temp_sharp_curr,Temp_sharp_consec);

% ---------------------- Tunable thresholds ----------------------
area_sum_threshold     = 0.1;   % Ratio threshold for multi-region comparison
maj_thresh             = 15;    % Major-axis threshold for large ellipses
del_var                = 0.05;  % Minimum overlap ratio for keeping candidates
allow_area_decr_thresh = 0.8;   % Lower bound when preferring smaller ellipse

% Decision flags
curr_flag   = 1;   % Replace ellipse in Comp_Bin_Im
consec_flag = 2;   % Keep existing ellipse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2A: Resolve competing regions

if ~isempty(Curr_Cell)
    parfor tt = 1:size(Curr_Cell,2)

        % Local copies for parfor safety
        Curr_Cell_TMP          = Curr_Cell;
        areas_sharp_curr_TMP   = areas_sharp_curr;
        areas_sharp_consec_TMP = areas_sharp_consec;

        % -------- ONE-TO-ONE MATCH --------
        if size(Curr_Cell_TMP{tt},2) == 1

            % Extract competing areas
            area_extract = [ ...
                areas_sharp_curr_TMP(Curr_Cell_TMP{tt}(1)), ...
                areas_sharp_consec_TMP(Curr_Cell_TMP{tt}(2)) ];

            % Delegate decision logic
            [choose_conf_idx{tt}] = UNIQUE_MATCH_ANALYSIS( ...
                area_extract,orig_comp,Pix_idx_consec,Curr_Cell_TMP, ...
                allow_area_decr_thresh,consec_flag,curr_flag, ...
                Clean_matrix,tt,Base_AREA_compare,Comp_Bin_Im);

        % -------- ONE-TO-MANY MATCH --------
        elseif size(Curr_Cell_TMP{tt},2) > 1

            % Handles potential oversegmentation cases
            [choose_conf_idx{tt}] = NONUNIQUE_MATCH_ANALYSIS( ...
                Curr_Cell_TMP,Pix_idx_consec,Pix_idx_curr,del_var, ...
                consec_flag,curr_flag,maj_thresh,area_sum_threshold,tt, ...
                Stats_Current,Stats_Consecutive, ...
                areas_sharp_curr_TMP,areas_sharp_consec_TMP);
        end
    end
end

% Re-declare flags for clarity in later logic
curr_flag   = 1;  % Replace ellipse
consec_flag = 2;  % Keep ellipse

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 3: Assemble corrected statistics structure

if ~isempty(Curr_Cell)
    choose_conf_idx_LMP = [choose_conf_idx{:}];

    % Selected regions from each slice
    MAIN_STATS_CURR = Stats_Current( ...
        unique(choose_conf_idx_LMP(2, ...
        ismember(choose_conf_idx_LMP(1,:),curr_flag))));

    MAIN_STATS_CONSEC = Stats_Consecutive( ...
        unique(choose_conf_idx_LMP(2, ...
        ismember(choose_conf_idx_LMP(1,:),consec_flag))));
else
    MAIN_STATS_CURR   = [];
    MAIN_STATS_CONSEC = [];
end

% Collect unmatched and ambiguous regions
MAIN_STATS_CURR_MISSING   = Stats_Current(Curr_uncaptured_cell);
MAIN_STATS_CURR_NONUNIQUE = Stats_Current(unsure_probe_vals);
MAIN_STATS_CONSEC_MISSING = Stats_Consecutive(Consec_uncaptured_cell);

% Final merged stats
STATS_CORRECTED_FINAL = [ ...
    MAIN_STATS_CURR; ...
    MAIN_STATS_CURR_MISSING; ...
    MAIN_STATS_CONSEC; ...
    MAIN_STATS_CONSEC_MISSING; ...
    MAIN_STATS_CURR_NONUNIQUE ];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 3B: Update Comp_Bin_Im and Base_AREA_compare
% Purpose:
%   Maintain a stable reference area set and avoid progressive shrinkage
%   of optimized segmentations.

if orig_comp == 1
    % -------- INITIAL COMPARISON --------
    if ~isempty(Curr_Cell)
        tot_current_chosen = choose_conf_idx_LMP(2, ...
            ismember(choose_conf_idx_LMP(1,:),curr_flag));
    else
        tot_current_chosen = [];
    end

    % Remaining regions
    tot_current_remain = [ ...
        find(~ismember(1:size(Stats_Current,1),tot_current_chosen)), ...
        Curr_uncaptured_cell];

    tot_consec_remain = Consec_uncaptured_cell;

    % Collect pixel indices
    Pix_idx_TOTAL = [ ...
        {Stats_Current(tot_current_remain).Pixel_IDX_List}, ...
        {Stats_Consecutive(tot_consec_remain).Pixel_IDX_List}];

    % Update reference areas
    Base_AREA_compare = [ ...
        areas_sharp_curr(tot_current_remain), ...
        areas_sharp_consec(tot_consec_remain)];

    % Populate label image
    for tt = 1:length(Pix_idx_TOTAL)
        Comp_Bin_Im(Pix_idx_TOTAL{tt}) = tt;
    end

elseif orig_comp == 0
    % -------- SUBSEQUENT COMPARISONS --------

    tot_current_remain = Curr_uncaptured_cell;
    tot_consec_remain  = Consec_uncaptured_cell;

    Pix_idx_TOTAL = [ ...
        {Stats_Current(tot_current_remain).Pixel_IDX_List}, ...
        {Stats_Consecutive(tot_consec_remain).Pixel_IDX_List}];

    Base_AREA_compare = [ ...
        Base_AREA_compare, ...
        areas_sharp_curr(tot_current_remain), ...
        areas_sharp_consec(tot_consec_remain)];

    % Append labels without overwriting existing ones
    if ~isempty(Pix_idx_TOTAL)
        for tt = 1:length(Pix_idx_TOTAL)
            if all(unique(Comp_Bin_Im(Pix_idx_TOTAL{tt})) == 0)
                Comp_Bin_Im(Pix_idx_TOTAL{tt}) = base_area_length + tt;
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 4: Final geometric correction

C_value = 0.90;
STATS_CORRECTED_FINAL = INTERSECTOR_CORRECTOR( ...
    STATS_CORRECTED_FINAL,Clean_matrix,C_value, ...
    columnsInImage,rowsInImage);

% ---------------------- OUTPUTS ----------------------
varargout{1}      = STATS_CORRECTED_FINAL;
varargout{end+1}  = Comp_Bin_Im;
varargout{end+1}  = Base_AREA_compare;

end




% =========================================================================
% UNIQUE_MATCH_ANALYSIS
% =========================================================================
% This function performs a one-to-one comparative analysis between ellipses
% from two segmentation results. It selects which ellipse to retain based
% on area comparison and greyscale intensity criteria.
%
% INPUTS (in order):
%   area_extract             - Array containing the areas of the two ellipses
%   orig_comp                - Flag indicating whether this is the first comparison
%                              between the first two images (1) or a subsequent
%                              comparison (0)
%   Pix_idx_consec_TMP       - Cell array of linear indices of ellipses for the
%                              consecutive image slice
%   Curr_Cell_TMP            - Cell array capturing indices of overlapping ellipses
%                              when one image is overlaid on another
%   allow_area_decr_thresh   - Allowable ratio of overlap between ellipses
%   consec_flag              - Flag indicating the previous ellipse should be kept
%   curr_flag                - Flag indicating the current ellipse should replace the previous
%   Clean_matrix_TMP         - Original greyscale image
%   tt                       - Current loop index
%   Base_AREA_compare_TMP    - Areas of the ellipses in the current optimized result
%   Comp_Bin_Im_TMP          - Current optimized segmentation result
%
% OUTPUTS:
%   choose_conf_idx          - Selected ellipse indices based on area and intensity criteria
% =========================================================================

function [varargout] = UNIQUE_MATCH_ANALYSIS(varargin)

% Assign inputs to readable variable names
area_extract             = varargin{1};
orig_comp                = varargin{2};
Pix_idx_consec_TMP       = varargin{3};
Curr_Cell_TMP            = varargin{4};
allow_area_decr_thresh   = varargin{5};
consec_flag              = varargin{6};
curr_flag                = varargin{7};
Clean_matrix_TMP         = varargin{8};
tt                       = varargin{9};
Base_AREA_compare_TMP    = varargin{10};
Comp_Bin_Im_TMP          = varargin{11};

%% ------------------------------------------------------------------------
% Step 1: Compare the areas of the two candidate ellipses
% -------------------------------------------------------------------------
if (area_extract(2) <= area_extract(1))
    
    if orig_comp == 1
        % Comparing first two images: use first area as reference
        Area_find = area_extract(1);
    else
        % Comparing a subsequent image to the optimized current image
        Area_find = Base_AREA_compare_TMP(nonzeros(Comp_Bin_Im_TMP(Pix_idx_consec_TMP{Curr_Cell_TMP{tt}(2)})));
    end
    
    % Determine which ellipse to choose based on allowable overlap ratio
    if (~isempty(Area_find)) && (length(unique(Area_find)) == 1) && ...
            ((area_extract(2) / unique(Area_find)) >= allow_area_decr_thresh)
        choose_conf_idx = [consec_flag ; Curr_Cell_TMP{tt}(2)];
    else
        choose_conf_idx = [curr_flag ; Curr_Cell_TMP{tt}(1)];
    end

%% ------------------------------------------------------------------------
% Step 2: Handle case where second area exceeds the first
% -------------------------------------------------------------------------
elseif area_extract(2) > area_extract(1)
    % Extract greyscale values of the second ellipse
    greyscale_consec = Clean_matrix_TMP(Pix_idx_consec_TMP{Curr_Cell_TMP{tt}(2)});
    
    set_thresh = 46250;   % Threshold intensity
    set_ratio  = 0.5;     % Fraction of pixels above threshold
    
    if (numel(find(greyscale_consec > set_thresh)) / numel(greyscale_consec)) >= set_ratio
        % Keep the consecutive ellipse if intensity criterion is met
        choose_conf_idx = [consec_flag ; Curr_Cell_TMP{tt}(2)];
    else
        % Otherwise, keep the current ellipse
        choose_conf_idx = [curr_flag ; Curr_Cell_TMP{tt}(1)];
    end
end         

%% ------------------------------------------------------------------------
% Step 3: Return the selected ellipse indices
% -------------------------------------------------------------------------
varargout{1} = choose_conf_idx;

end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [varargout] = NONUNIQUE_MATCH_ANALYSIS(varargin)
%--------------------------------------------------------------------------
% NONUNIQUE_MATCH_ANALYSIS
%
% PURPOSE:
%   Resolves conflicts in ellipse matching between two consecutive slices
%   when the relationship is not one-to-one (i.e., one-to-many or many-to-one).
%   The function decides whether to KEEP ellipses from the consecutive slice
%   or REPLACE them with ellipses from the current slice based on overlap,
%   area ratios, and major axis thresholds.
%
% MATCHING CASES HANDLED:
%   • 2:1  (many current → one consecutive)
%   • 1:2  (one current → many consecutive)
%   • many:many (no clear dominance)
%--------------------------------------------------------------------------

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% VARIABLE DEFINITIONS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Curr_Cell_TMP            -> Cell array storing matched ellipse indices
% Pix_idx_consec_TMP       -> Linear pixel indices (consecutive slice)
% Pix_idx_curr_TMP         -> Linear pixel indices (current slice)
% del_var                  -> Overlap deletion threshold
% consec_flag              -> Flag indicating KEEP from consecutive slice
% curr_flag                -> Flag indicating REPLACE with current slice
% maj_thresh               -> Major axis length threshold
% area_sum_threshold       -> Area ratio threshold for decision
% tt                       -> Loop index into Curr_Cell_TMP
% Stats_Current_TMP        -> Regionprops metadata (current slice)
% Stats_Consecutive_TMP    -> Regionprops metadata (consecutive slice)
% areas_sharp_curr_TMP     -> Ellipse areas (current slice)
% areas_sharp_consec_TMP   -> Ellipse areas (consecutive slice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-----------------------%
% Unpack input arguments
%-----------------------%
Curr_Cell_TMP              = varargin{1};
Pix_idx_consec_TMP         = varargin{2};
Pix_idx_curr_TMP           = varargin{3};
del_var                    = varargin{4};
consec_flag                = varargin{5};
curr_flag                  = varargin{6};
maj_thresh                 = varargin{7};
area_sum_threshold         = varargin{8};
tt                         = varargin{9};
Stats_Current_TMP          = varargin{10};
Stats_Consecutive_TMP      = varargin{11};
areas_sharp_curr_TMP       = varargin{12};
areas_sharp_consec_TMP     = varargin{13};

%==========================================================================
% CASE 1: MANY (current) → ONE (consecutive)  [2:1 matching]
%==========================================================================
if  numel(unique(Curr_Cell_TMP{tt}(1,:))) ~= 1 && ...
    numel(unique(Curr_Cell_TMP{tt}(2,:))) == 1
           
    % Extract pixel index sets
    unique_fib_idx     = Pix_idx_consec_TMP{unique(Curr_Cell_TMP{tt}(2,:))};
    non_unique_fib_idx = Pix_idx_curr_TMP(Curr_Cell_TMP{tt}(1,:));
           
    % Identify index groups
    nonunique_indices  = unique(Curr_Cell_TMP{tt}(1,:));
    unique_indices     = unique(Curr_Cell_TMP{tt}(2,:));
           
    % Compute total areas for comparison
    area_nonunique     = sum(areas_sharp_curr_TMP(nonunique_indices));
    area_unique        = areas_sharp_consec_TMP(unique_indices);
           
    % Remove weak-overlap matches
    del_tmp = cell2mat(cellfun(@(x,y) ...
        numel(find(ismember(x,unique_fib_idx)))/numel(x), ...
        non_unique_fib_idx,'uni',0)) < del_var;
           
    if ~isempty(nonzeros(del_tmp))
        Curr_Cell_TMP{tt}(:,del_tmp) = [];
    end
           
    % Decision based on major axis and area ratio
    if Stats_Consecutive_TMP(unique_indices).MajorAxisLength > maj_thresh  
        if ((area_nonunique/area_unique) < area_sum_threshold)
            % Prefer consecutive ellipse (KEEP)
            choose_conf_idx = [consec_flag*ones(1,size(Curr_Cell_TMP{tt},2)); ...
                               Curr_Cell_TMP{tt}(consec_flag,:)];
        else
            % Prefer current ellipses (REPLACE)
            choose_conf_idx = [curr_flag*ones(1,size(Curr_Cell_TMP{tt},2)); ...
                               Curr_Cell_TMP{tt}(curr_flag,:)];
        end
    elseif Stats_Consecutive_TMP(unique_indices).MajorAxisLength < maj_thresh  
        % Consecutive ellipse is too small → KEEP consecutive
        choose_conf_idx = [consec_flag*ones(1,size(Curr_Cell_TMP{tt},2)); ...
                           Curr_Cell_TMP{tt}(consec_flag,:)];
    end
           
%==========================================================================
% CASE 2: ONE (current) → MANY (consecutive)  [1:2 matching]
%==========================================================================
elseif numel(unique(Curr_Cell_TMP{tt}(1,:))) == 1 && ...
       numel(unique(Curr_Cell_TMP{tt}(2,:))) ~= 1
           
    % Extract pixel index sets
    unique_fib_idx       = Pix_idx_curr_TMP{unique(Curr_Cell_TMP{tt}(1,:))};
    non_unique_fib_idx   = Pix_idx_consec_TMP(Curr_Cell_TMP{tt}(2,:));
           
    % Identify index groups
    unique_indices       = unique(Curr_Cell_TMP{tt}(1,:));
    nonunique_indices    = unique(Curr_Cell_TMP{tt}(2,:));
         
    % Compute areas
    area_nonunique       = sum(areas_sharp_consec_TMP(nonunique_indices));
    area_unique          = areas_sharp_curr_TMP(unique_indices);
           
    % Remove weak-overlap matches
    del_tmp = cell2mat(cellfun(@(x,y) ...
        numel(find(ismember(x,unique_fib_idx)))/numel(x), ...
        non_unique_fib_idx,'uni',0)) < del_var;   

    if ~isempty(nonzeros(del_tmp))
        Curr_Cell_TMP{tt}(:,del_tmp) = [];
    end  
           
    % Default preference → consecutive
    choose_conf_idx = [consec_flag*ones(1,size(Curr_Cell_TMP{tt},2)); ...
                       Curr_Cell_TMP{tt}(consec_flag,:)];
           
    % Refine decision using geometry
    if Stats_Current_TMP(unique_indices).MajorAxisLength > maj_thresh
        if ((area_nonunique/area_unique) < area_sum_threshold)
            % Prefer current ellipse
            choose_conf_idx = [curr_flag*ones(1,size(Curr_Cell_TMP{tt},2)); ...
                               Curr_Cell_TMP{tt}(curr_flag,:)];
        else
            % Prefer consecutive ellipses
            choose_conf_idx = [consec_flag*ones(1,size(Curr_Cell_TMP{tt},2)); ...
                               Curr_Cell_TMP{tt}(consec_flag,:)];
        end
               
    elseif Stats_Current_TMP(unique_indices).MajorAxisLength < maj_thresh
        % Current ellipse too small → choose current
        choose_conf_idx = [curr_flag*ones(1,size(Curr_Cell_TMP{tt},2)); ...
                           Curr_Cell_TMP{tt}(curr_flag,:)];
    end
 
%==========================================================================
% CASE 3: MANY ↔ MANY (ambiguous matching)
%==========================================================================
elseif numel(unique(Curr_Cell_TMP{tt}(1,:))) ~= 1 && ...
       numel(unique(Curr_Cell_TMP{tt}(2,:))) ~= 1
           
    % No strong dominance — keep current indices
    choose_conf_idx = [curr_flag*ones(1,size(unique(Curr_Cell_TMP{tt}(1,:)),2)); ...
                       unique(Curr_Cell_TMP{tt}(1,:))];
end

%--------------------------------------------------------------------------
% Output
%--------------------------------------------------------------------------
varargout{1} = choose_conf_idx;            
end