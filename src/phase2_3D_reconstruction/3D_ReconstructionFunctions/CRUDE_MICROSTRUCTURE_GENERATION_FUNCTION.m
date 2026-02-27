%% =============================================================================
% Current Update: 3/16/2018
% =============================================================================
% CRUDE_MICROSTRUCTURE_GENERATION_FUNCTION
%
% Description:
% This function generates a preliminary 3D microstructure by stacking
% ellipses fitted to segmented regions of 2D slices. It performs:
%   1. Removal of redundant metadata
%   2. Similar centroid comparison for initial reconstruction
%   3. Tracking of ellipse indexes across slices
%   4. Extraction of constituent properties
%   5. Creation of pixel lists for each fiber
%   6. Conversion to linear indices for crude 3D fiber volumes
%
% INPUTS (varargin):
% ------------------
% SLICE_REGIONS_STATS_STRUCTURE : Structure array with metadata of fitted ellipses per slice
% Centroid_cell                 : Cell array of centroids for each slice
% numberOfImageFiles            : Number of slices in the dataset
% size_length                   : Size of the 2D greyscale image slices
% iter_limit                     : Maximum slice distance to stack ellipses
%
% OUTPUTS (varargout):
% --------------------
% 1. query                    : String 'done' as a process flag
% 2. properties               : Metadata of constituent ellipses
% 3. index_crude_fibers       : Linear indices of each fiber in 3D volume
% 4. cat_size_const_ellipses  : Number of pixels in each ellipse
%% =============================================================================

function [varargout]= CRUDE_MICROSTRUCTURE_GENERATION_FUNCTION(varargin)

% --------------------------
% Extract input arguments
% --------------------------
SLICE_REGIONS_STATS_STRUCTURE = varargin{1};
Centroid_cell                 = varargin{2};
numberOfImageFiles            = varargin{3};
size_length                   = varargin{4};
iter_limit                    = varargin{5};

%% ============================================================================
% STEP 2: Remove repeated metadata to ensure uniqueness
% ============================================================================
SLICE_UPDATE    = SLICE_REGIONS_STATS_STRUCTURE;
CENTROID_UPDATE = Centroid_cell;

% rept_correction_TT removes redundant entries in the slice library
[SLICE_UPDATE, CENTROID_UPDATE] = rept_correction_TT(SLICE_UPDATE, CENTROID_UPDATE, size_length);

%% ============================================================================
% STEP 3: Similar centroid comparison
% ============================================================================
% Compares centroids across slices to identify potentially connected ellipses
clearvars sieve_parameter_unique_CELL remnant_ellipses track_MAIN_CELL properties cent_info_diff_check

disp('BEGIN SIMIL_CENTROID_COMPARISON')
[sieve_parameter_unique_CELL, remnant_ellipses] = SIMIL_CENTROID_COMPARISON_UPDATED_TT(numberOfImageFiles, SLICE_UPDATE, CENTROID_UPDATE);
disp('END SIMIL_CENTROID_COMPARISON')

%% ============================================================================
% STEP 4: Index tracking
% ============================================================================
% Tracks the indices of ellipses across slices to facilitate stacking
disp('BEGIN INDEX_TRACK_ESTIMATOR')
numtrack = size(SLICE_UPDATE, 2) - 1; 
[track_MAIN_CELL, ~] = INDEX_TRACK_ESTIMATOR(numberOfImageFiles, remnant_ellipses, sieve_parameter_unique_CELL, iter_limit, numtrack);
disp('END INDEX_TRACK_ESTIMATOR')

%% ============================================================================
% STEP 5: Extract properties for each constituent ellipse
% ============================================================================
% Computes geometric and positional properties for each fiber
disp('BEGIN FIBER_CONSTITUENT_PROPERTIES_ESTIMATOR')
[properties] = FIBER_CONSTITUENT_PROPERTIES_ESTIMATOR(track_MAIN_CELL, SLICE_UPDATE);
disp('END FIBER_CONSTITUENT_PROPERTIES_ESTIMATOR')

%% ============================================================================
% STEP 6: Generate pixel lists for each ellipse
% ============================================================================
% Each cell contains [x y slice] coordinates of all pixels in an ellipse

clearvars pix_list_cell
size_const_ellipses = cell(1, size(properties,2));
pix_list_cell       = cell(1, size(properties,2));

for ii = 1:size(properties,2)
    if ~isempty(properties{ii})
        % Count number of pixels per ellipse
        size_const_ellipses{ii} = cellfun(@(x)size(x,1), properties{ii}, 'uni', 0);
        
        % Construct pixel lists with slice index
        for kk = 1:size(properties{ii},2)
            if ii == 1
                % Slice index starts at 1 for first group
                pix_list_cell{ii}{kk} = cellfun(@(x,y)[x y*ones(size(x,1),1)], ...
                                                {properties{ii}{kk}.Pixel_List}, ...
                                                num2cell(1:length({properties{ii}{kk}.Pixel_List})), ...
                                                'uni',0);
            else 
                % Subsequent slices get adjusted slice indices
                pix_list_cell{ii}{kk} = cellfun(@(x,y)[x y*ones(size(x,1),1)], ...
                                                {properties{ii}{kk}.Pixel_List}, ...
                                                num2cell(ii:ii+length({properties{ii}{kk}.Pixel_List})-1), ...
                                                'uni',0);
            end
        end
    end
end

% Concatenate all pixel lists
cat_pix_list_cell = cell(1,size(properties,2));
for ii = 1:size(properties,2)   
    for kk = 1:size(properties{ii},2)
        cat_pix_list_cell{ii}{kk} = vertcat(pix_list_cell{ii}{kk}{:});
    end    
end

MAIN_cat_pix_list_cell  = horzcat(cat_pix_list_cell{:});
cat_size_const_ellipses = cell2mat(horzcat(size_const_ellipses{:}));

%% ============================================================================
% STEP 7: Convert pixel coordinates to linear indices in the 3D volume
% ============================================================================
% Allows direct indexing of fibers in the 3D volume for reconstruction
index_crude_fibers = cell(1, size(MAIN_cat_pix_list_cell,2));

for pp = 1:size(MAIN_cat_pix_list_cell,2)
    index_crude_fibers{pp} = sub2ind(size_length, ...
                                     MAIN_cat_pix_list_cell{pp}(:,2), ...
                                     MAIN_cat_pix_list_cell{pp}(:,1), ...
                                     MAIN_cat_pix_list_cell{pp}(:,3));   
end

% --------------------------
% Output
% --------------------------
query = 'done'; % Flag indicating completion

varargout{1} = query;
varargout{end+1} = properties;
varargout{end+1} = index_crude_fibers;
varargout{end+1} = cat_size_const_ellipses;