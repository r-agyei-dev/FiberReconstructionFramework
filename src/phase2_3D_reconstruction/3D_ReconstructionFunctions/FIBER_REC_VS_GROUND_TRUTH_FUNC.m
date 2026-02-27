%% =============================================================================
% FIBER_REC_VS_GROUND_TRUTH_FUNC
% =============================================================================
% Description:
% This function prepares labeled volumetric images of reconstructed fibers
% and their corresponding slice indices, then applies the EPAD stitching
% algorithm to resolve overlapping fibers and identify contiguous segments.
%
% INPUTS (via varargin):
% ----------------------
% Linear_Index_Final   : Cell array of linear indices for each fiber
% height_variable_MAIN : Slice number associated with each fiber pixel
% Pixels_NON_BIN       : Logical mask indicating non-connected pixels
% size_length          : Size of the 3D image volume
% save_path_extra_file : (unused here, legacy for saving outputs)
% load_variable        : (unused)
% save_path            : (unused)
% new_idx              : (unused)
%
% OUTPUTS:
% --------
% varargout{1} : split_array_slice       (slice groupings per fiber)
% varargout{2} : split_array_consec_int  (contiguous index groupings)
%
% NOTES:
% ------
% - Builds two labeled volumes:
%     IMAGE_VOL_1 → fiber IDs
%     IMAGE_VOL_2 → slice numbers
% - Uses EPAD algorithm to resolve overlaps
% - Memory clearing is intentionally preserved
%% =============================================================================

function [varargout] = FIBER_REC_VS_GROUND_TRUTH_FUNC(varargin)

%% ============================================================================
% STEP 0: Parse inputs
% ============================================================================
Linear_Index_Final   = varargin{1};
height_variable_MAIN = varargin{2};
Pixels_NON_BIN       = varargin{3};
size_length          = varargin{4};
save_path_extra_file = varargin{5}; %#ok<NASGU>
load_variable        = varargin{6}; %#ok<NASGU>
save_path            = varargin{7}; %#ok<NASGU>
new_idx              = varargin{8}; %#ok<NASGU>

tic

%% ============================================================================
% STEP 1: Instantiate labeled image volumes
% ============================================================================
disp('Begin Image instantiations')

% Volume storing fiber IDs
IMAGE_VOL_1 = zeros(size_length);

% Volume storing slice numbers
IMAGE_VOL_2 = zeros(size_length);

% Preallocate helper cell arrays
curr_values = cell(1, length(Linear_Index_Final));
succ_values = cell(1, length(Linear_Index_Final));

%% ============================================================================
% STEP 2: Populate volumes with fiber information
% ============================================================================
for ii = 1:length(Linear_Index_Final)

    % Label fibers in IMAGE_VOL_1 using fiber index
    IMAGE_VOL_1(Linear_Index_Final{ii}) = ii;
    curr_values{ii} = IMAGE_VOL_1(Linear_Index_Final{ii});

    % Label slice numbers in IMAGE_VOL_2
    IMAGE_VOL_2(Linear_Index_Final{ii}) = height_variable_MAIN{ii};

    % Remove non-connected regions by setting to zero
    IMAGE_VOL_2(Linear_Index_Final{ii}(Pixels_NON_BIN(Linear_Index_Final{ii}))) = 0;

    % Store successor values
    succ_values{ii} = IMAGE_VOL_2(Linear_Index_Final{ii});

end

disp('End Image instantiations')

%% ============================================================================
% STEP 3: Clear large variables to relieve memory
% ============================================================================
clearvars height_variable_MAIN Linear_Index_Final Pixels_NON_BIN
clearvars Pixels_NON_BIN   % (intentional duplicate preserved)

%% ============================================================================
% STEP 4: Run EPAD stitching algorithm
% ============================================================================
disp('Use the EPAD Algorithm and solve')

% Compute exponent used for unique encoding
exponent = 10^(numel(num2str(max(IMAGE_VOL_2(:)))) + 2);

% Run EPAD fiber stitching
[curr_succ_multiple, split_array_consec_int, split_array_slice] = ...
    EPAD_FUNCTION_CODE_Fiber_Stitch_ver(curr_values, succ_values, exponent);

% Clear large image volumes
clearvars IMAGE_VOL_1 IMAGE_VOL_2

toc

%% ============================================================================
% STEP 5: Assign outputs
% ============================================================================
varargout{1}     = split_array_slice;
varargout{end+1} = split_array_consec_int;

end