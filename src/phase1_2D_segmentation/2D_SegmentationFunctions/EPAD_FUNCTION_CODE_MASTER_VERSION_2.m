%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version created on: 9/2/2018
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%--------------------------------------------------------------------------
% EPAD_FUNCTION_CODE_MASTER_VERSION_2
%
% PURPOSE:
%   Performs fast overlap analysis between two labeled 2D volumes
%   (current vs successive). The EPAD arithmetic encoding is used to
%   efficiently detect intersecting region indices.
%
% OUTPUTS:
%   1) Curr_Cell                -> cell array of matched curr/succ indices
%   2) Curr_uncaptured_cell     -> curr indices with no successor overlap
%   3) Consec_uncaptured_cell   -> successor indices with no curr overlap
%   4) unsure_probe_vals        -> ambiguous matches removed from analysis
%
% METHOD (EPAD):
%   Overlaps are detected using arithmetic encoding based on:
%     • Exponentiation
%     • Premultiplication
%     • Addition
%     • Division
%--------------------------------------------------------------------------

function [varargout] = EPAD_FUNCTION_CODE_MASTER_VERSION_2(varargin)

%--------------------------------------------------------------------------
% Unpack inputs (labeled images)
%--------------------------------------------------------------------------
curr_vol = varargin{1};   % current slice label map
succ_vol = varargin{2};   % successive slice label map

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Encode overlaps (current → successive)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pre_vol  = curr_vol;
post_vol = succ_vol;

% Create scaling factor large enough to separate labels numerically
exponent = 10^(numel(num2str(max(post_vol(:)))) + 2);

% Arithmetic encoding of pairwise overlaps
pre_vol_pre  = pre_vol * exponent;
curr_succ_SUM = pre_vol_pre + post_vol;

% Count occurrences of each encoded overlap value
clearvars p_idx_mod p_idx_mod_1 curr_succ_overlap
[p_idx_mod(:,2),p_idx_mod(:,1)] = ...
    hist(nonzeros(curr_succ_SUM(:)), unique(nonzeros(curr_succ_SUM(:))));

% Remove values that do not represent valid overlaps
p_idx_mod(p_idx_mod(:,1) < exponent,:) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Decode overlap table
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Column 1: current index
curr_succ_overlap(:,1) = floor(p_idx_mod(:,1) / exponent);

% Column 2: successive index
curr_succ_overlap(:,2) = mod(p_idx_mod(:,1), exponent);

% Capture current fibers intersecting only background
potential_curr_empty = curr_succ_overlap(curr_succ_overlap(:,2) == 0, 1);

% Remove background intersections from overlap table
curr_succ_overlap(curr_succ_overlap(:,2) == 0,:) = [];

%% Group matches by current and successor indices
table_fin = curr_succ_overlap(:,1:2);

%--------------------------------------------------------------------------
% Build grouping from CURRENT perspective
%--------------------------------------------------------------------------
unique_fib_idx = nonzeros(unique(table_fin(:,1)));
MAIN_COUNT_probe12_TMP_1 = cell(1,numel(unique_fib_idx));

for ii = 1:numel(unique_fib_idx)
    MAIN_COUNT_probe12_TMP_1{ii} = ...
        (table_fin(ismember(table_fin(:,1), unique_fib_idx(ii)),:))';
end

%--------------------------------------------------------------------------
% Build grouping from SUCCESSOR perspective
%--------------------------------------------------------------------------
unique_fib_idx = nonzeros(unique(table_fin(:,2)));
MAIN_COUNT_query12_TMP_1 = cell(1,numel(unique_fib_idx));

for ii = 1:numel(unique_fib_idx)
    MAIN_COUNT_query12_TMP_1{ii} = ...
        (table_fin(ismember(table_fin(:,2), unique_fib_idx(ii)),:))';
end

% Remove trivial one-to-one groups from successor side
MAIN_COUNT_query12_TMP_1( ...
    cell2mat(cellfun(@(x)size(x,2)==1, ...
    MAIN_COUNT_query12_TMP_1,'uni',0))) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove probe groups that are singletons but belong to larger successor groups
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~isempty(MAIN_COUNT_probe12_TMP_1) && ~isempty(MAIN_COUNT_query12_TMP_1)

    query_out = [MAIN_COUNT_query12_TMP_1{:}];

    for ii = 1:size(MAIN_COUNT_probe12_TMP_1,2)
        if size(MAIN_COUNT_probe12_TMP_1{ii},2) == 1 && ...
           ismember(MAIN_COUNT_probe12_TMP_1{ii}(2), unique(query_out(2,:)))
            MAIN_COUNT_probe12_TMP_1{ii} = [];
        end
    end

    % Remove emptied cells
    MAIN_COUNT_probe12_TMP_1( ...
        cell2mat(cellfun(@(x)isempty(x), ...
        MAIN_COUNT_probe12_TMP_1,'uni',0))) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Repeat encoding in reverse (successive → current)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pre_vol  = succ_vol;
post_vol = curr_vol;

exponent = 10^(numel(num2str(max(post_vol(:)))) + 2);

pre_vol_pre   = pre_vol * exponent;
succ_curr_SUM = pre_vol_pre + post_vol;

clearvars p_idx_mod p_idx_mod_1 curr_succ_overlap
[p_idx_mod(:,2),p_idx_mod(:,1)] = ...
    hist(nonzeros(succ_curr_SUM(:)), unique(nonzeros(succ_curr_SUM(:))));

p_idx_mod(p_idx_mod(:,1) < exponent,:) = [];

% Decode reverse overlaps
succ_curr_overlap(:,1) = floor(p_idx_mod(:,1) / exponent);
succ_curr_overlap(:,2) = mod(p_idx_mod(:,1), exponent);

% Successor fibers intersecting only background
potential_succ_empty = succ_curr_overlap(succ_curr_overlap(:,2) == 0, 1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Remove ambiguous many-to-many conflicts
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Curr_Cell = [MAIN_COUNT_probe12_TMP_1 MAIN_COUNT_query12_TMP_1];

if ~isempty(MAIN_COUNT_probe12_TMP_1) && ~isempty(MAIN_COUNT_query12_TMP_1)

    probe_lump = [MAIN_COUNT_probe12_TMP_1{:}];
    query_lump = [MAIN_COUNT_query12_TMP_1{:}];

    % Remove already-accounted empty indices
    potential_curr_empty(ismember(potential_curr_empty, ...
        [probe_lump(1,:) query_lump(1,:)])) = [];

    potential_succ_empty(ismember(potential_succ_empty, ...
        [query_lump(2,:) probe_lump(2,:)])) = [];

    % Identify ambiguous (non-unique) mappings
    non_unique_probe_1 = ...
        probe_lump(1, ismember(probe_lump(1,:), query_lump(1,:)));

    non_unique_query = ...
        query_lump(2, ismember(query_lump(2,:), probe_lump(2,:)));

    non_unique_probe_2 = ...
        query_lump(1, ismember(query_lump(2,:), non_unique_query));

    unsure_probe_vals = unique([non_unique_probe_1 non_unique_probe_2]);

    % Remove ambiguous groups from final matching
    Curr_Cell( ...
        cell2mat(cellfun(@(x)any(ismember(x(1,:),unsure_probe_vals)), ...
        Curr_Cell,'uni',0))) = [];

else
    unsure_probe_vals = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Final outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Curr_uncaptured_cell   = potential_curr_empty';
Consec_uncaptured_cell = potential_succ_empty';

varargout{1}     = Curr_Cell;
varargout{end+1} = Curr_uncaptured_cell;
varargout{end+1} = Consec_uncaptured_cell;
varargout{end+1} = unsure_probe_vals;

end