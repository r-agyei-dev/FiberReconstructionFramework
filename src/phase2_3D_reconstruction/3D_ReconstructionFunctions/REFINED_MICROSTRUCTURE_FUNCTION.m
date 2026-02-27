% =========================================================================
% PURPOSE:
% Perform refined 3D microstructure reconstruction and labeling of fibers.
%
% This function:
%   • Builds labeled 3D fiber volumes
%   • Computes fiber height indices
%   • Extracts centroids and out-of-plane angles
%   • Saves key metadata for later analysis
%
% The processing path depends on the orientation specified by load_variable.
%
% OUTPUTS:
%   varargout{1} → status string
%   varargout{2} → missing fiber indices
%   varargout{3} → reconstructed labeled volume
% =========================================================================

function [varargout] = REFINED_MICROSTRUCTURE_FUNCTION(varargin)

% -------------------------------------------------------------------------
% INPUT ARGUMENTS
% -------------------------------------------------------------------------
Linear_Index_Final       = varargin{1};   % Fiber voxel linear indices
orig_size_length         = varargin{2};   % Original tomogram size
size_length              = varargin{3};   % Working volume size
Pixel_Idx_Per_Slice      = varargin{4};   % Pixel indices per slice
load_variable            = varargin{5};   % Orientation selector (1,2,3)
FINAL_CORRECTED_CENTROID = varargin{6};   % Precomputed centroids
tomo_dataset_idx         = varargin{7};   % Dataset identifier
save_path                = varargin{8};   % Output directory
height_variable_MAIN     = varargin{9};   % Slice height information
crop_region              = varargin{10};  % Cropping limits

%% ------------------------------------------------------------------------
% NOTE:
% Fibers with high out-of-plane angle may be filtered later.
% -------------------------------------------------------------------------

tic

% -------------------------------------------------------------------------
% Processing configuration flags
% -------------------------------------------------------------------------
hdfsave         = 1;    % Enable HDF saving inside called function
sieve_indicator = 1;    % Enable sieving/filtering of fibers
angular_thresh  = 60;   % Out-of-plane angle threshold (degrees)

% -------------------------------------------------------------------------
% Output file name prefixes
% -------------------------------------------------------------------------
Temp_vol_prefix              = 'REFINED_FIBERS_UNSIEVEDn';
Idx_Height_prefix            = 'Idx_Height_UNSIEVED';
Center_Idx_prefix            = 'Center_Idx_UNSIEVED';
Slice_Range_prefix           = 'Slice_Range_UNSIEVED';
Centroid_Coordiantes_prefix  = 'Cent_Coord_UNSIEVED';
FIB_PIX_CORD_prefix          = 'FIB_PIX_CORD_UNSIEVED';
Angle_prefix                 = 'Out_Angle_UNSEIVED';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====================== ORIENTATION-SPECIFIC PROCESSING ==================
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% -------------------------------------------------------------------------
% CASE 1: XYZ orientation
% -------------------------------------------------------------------------
if load_variable == 1

    % Build labeled volume and associated metadata
    [Temp_volume_xyz_v2_UP,Idx_Height_xyz_v2_UP,...
     check_mat_center_GLOBAL_VEC_xyz_v2_UP,...
     zh_xyz,FINAL_CORRECTED_CENTROID_xyz,...
     out_of_plane_angle_xyz,missing_idx_xyz,...
     FIB_PIX_CORD_xyz] = ...
        THREE_D_LABELLER_HDF_CREATOR_DIRECT_UODATE_v2_TT( ...
            Linear_Index_Final,orig_size_length,size_length,hdfsave,...
            Pixel_Idx_Per_Slice,load_variable,FINAL_CORRECTED_CENTROID,...
            tomo_dataset_idx,save_path,sieve_indicator,angular_thresh,...
            height_variable_MAIN,crop_region);

    % Store outputs in common variables
    TEMP_MAT   = Temp_volume_xyz_v2_UP;
    missing_idx = missing_idx_xyz;

    % Save key metadata
    save(fullfile(save_path,...
        sprintf('%s%s',Center_Idx_prefix,'_xyz.mat')),...
        'check_mat_center_GLOBAL_VEC_xyz_v2_UP','-v7.3');

    save(fullfile(save_path,...
        sprintf('%s%s',Centroid_Coordiantes_prefix,'_xyz.mat')),...
        'FINAL_CORRECTED_CENTROID_xyz','-v7.3');

% -------------------------------------------------------------------------
% CASE 2: ZXY orientation
% -------------------------------------------------------------------------
elseif load_variable == 2

    [Temp_volume_zxy_v2_UP,Idx_Height_zxy_v2_UP,...
     check_mat_center_GLOBAL_VEC_zxy_v2_UP,...
     zh_zxy,FINAL_CORRECTED_CENTROID_zxy,...
     out_of_plane_angle_zxy,missing_idx_zxy,...
     FIB_PIX_CORD_zxy] = ...
        THREE_D_LABELLER_HDF_CREATOR_DIRECT_UODATE_v2_TT( ...
            Linear_Index_Final,orig_size_length,size_length,hdfsave,...
            Pixel_Idx_Per_Slice,load_variable,FINAL_CORRECTED_CENTROID,...
            tomo_dataset_idx,save_path,sieve_indicator,angular_thresh,...
            height_variable_MAIN,crop_region);

    TEMP_MAT   = Temp_volume_zxy_v2_UP;
    missing_idx = missing_idx_zxy;

    save(fullfile(save_path,...
        sprintf('%s%s',Center_Idx_prefix,'_zxy.mat')),...
        'check_mat_center_GLOBAL_VEC_zxy_v2_UP','-v7.3');

    save(fullfile(save_path,...
        sprintf('%s%s',Centroid_Coordiantes_prefix,'_zxy.mat')),...
        'FINAL_CORRECTED_CENTROID_zxy','-v7.3');

% -------------------------------------------------------------------------
% CASE 3: YZX orientation
% -------------------------------------------------------------------------
elseif load_variable == 3

    tic
    [Temp_volume_yzx_v2_UP,Idx_Height_yzx_v2_UP,...
     check_mat_center_GLOBAL_VEC_yzx_v2_UP,...
     zh_yzx,FINAL_CORRECTED_CENTROID_yzx,...
     out_of_plane_angle_yzx,missing_idx_yzx,...
     FIB_PIX_CORD_yzx] = ...
        THREE_D_LABELLER_HDF_CREATOR_DIRECT_UODATE_v2_TT( ...
            Linear_Index_Final,orig_size_length,size_length,hdfsave,...
            Pixel_Idx_Per_Slice,load_variable,FINAL_CORRECTED_CENTROID,...
            tomo_dataset_idx,save_path,sieve_indicator,angular_thresh,...
            height_variable_MAIN,crop_region);
    toc

    TEMP_MAT   = Temp_volume_yzx_v2_UP;
    missing_idx = missing_idx_yzx;

    save(fullfile(save_path,...
        sprintf('%s%s',Center_Idx_prefix,'_yzx.mat')),...
        'check_mat_center_GLOBAL_VEC_yzx_v2_UP','-v7.3');

    save(fullfile(save_path,...
        sprintf('%s%s',Centroid_Coordiantes_prefix,'_yzx.mat')),...
        'FINAL_CORRECTED_CENTROID_yzx','-v7.3');

end

toc

% -------------------------------------------------------------------------
% OUTPUT ASSIGNMENT
% -------------------------------------------------------------------------
stat_report = 'QUICK_Property_Function_C';

varargout{1} = stat_report;   % Status message
varargout{2} = missing_idx;   % Missing fiber indices
varargout{3} = TEMP_MAT;      % Final labeled volume