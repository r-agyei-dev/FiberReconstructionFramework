% ST_CHOOSE_OVER_OR_QUERY
% This function selects the best fiber centroid to use when multiple candidate
% fibers are present. It considers both the "query" fiber and the "overlapping" fiber,
% evaluating direction vectors, angular thresholds, and centroid positions.  
% The function handles four main cases:
%   1. Both query_var_over_segTMP and query_var are non-empty
%   2. query_var is non-empty and query_var_over_segTMP is empty
%   3. query_var is empty and query_var_over_segTMP is non-empty
%   4. Both query_var and query_var_over_segTMP are empty
%
% Inputs:
%   query_var_over_segTMP      - Candidate fiber index from overlapping segmentation
%   query_var                  - Candidate fiber index from query segmentation
%   direction_vectors           - Direction vectors of fibers
%   d_vec_angles_overTMP        - Precomputed angular differences between fibers
%   Centroid_MID_cell_LUMP_ORIG_UP - Original centroids for fibers
%   Centroid_MID_cell_LUMP_TMP - Temporary centroids used for comparison
%   z_TMP                       - Cell containing z-positions for fibers
%   z                           - Updated z-positions for fibers
%   c_f                         - Current fiber index
%   STATS_Pix_LOC_LUMP_UP       - Pixel statistics for fibers
%   ang_thresh                  - Angular threshold for matching fibers
%   main_res                     - Main results containing fiber coordinates
%
% Outputs:
%   NEW_CENTROID               - Selected or updated centroid
%   c_f                        - Updated current fiber index
%   z                          - Updated z positions
%   z_TMP                      - Updated temporary z positions
%   Centroid_MID_cell_LUMP_ORIG_UP - Updated original centroid cell array
%   Centroid_MID_cell_LUMP_TMP - Updated temporary centroid cell array
%   STATS_Pix_LOC_LUMP_UP      - Updated pixel statistics for fibers

function [varargout]   =   ST_choose_over_or_query(varargin)

query_var_over_segTMP                 = varargin{1};
query_var                             = varargin{2};
direction_vectors                     = varargin{3};
d_vec_angles_overTMP                  = varargin{4};
Centroid_MID_cell_LUMP_ORIG_UP        = varargin{5};
Centroid_MID_cell_LUMP_TMP            = varargin{6};
z_TMP                                 = varargin{7};
z                                     = varargin{8};
c_f                                   = varargin{9};
STATS_Pix_LOC_LUMP_UP                 = varargin{10};
ang_thresh                            = varargin{11};
main_res                              = varargin{12};


% CASE 1: Both query and overlapping fiber candidates exist
if ~isempty(query_var_over_segTMP) &&  ~isempty(query_var)
    
    % Compute angular difference and select the best candidate
    [~,final_check_tmp] =  min([acosd(dot(direction_vectors{c_f},direction_vectors{query_var})/(norm(direction_vectors{c_f})*norm(direction_vectors{query_var})))   d_vec_angles_overTMP]);
    final_check_idx = final_check_tmp(1);
    
    if final_check_idx == 1
        % Use query fiber
        [NEW_CENTROID] = probe_var_correction_up(Centroid_MID_cell_LUMP_ORIG_UP,c_f,query_var,main_res,z_TMP);
        c_f = query_var;
    else
        % Use overlapping fiber and update all relevant variables
        [NEW_CENTROID,Centroid_MID_cell_LUMP_ORIG_UP{query_var_over_segTMP},...
            z{query_var_over_segTMP},STATS_Pix_LOC_LUMP_UP{query_var_over_segTMP}] = Sticth_Shared_fibers(Centroid_MID_cell_LUMP_ORIG_UP,z,c_f,...
                                                                                                         query_var_over_segTMP,STATS_Pix_LOC_LUMP_UP);
        c_f = query_var_over_segTMP;
        z_TMP{c_f} = z{query_var_over_segTMP};       % Update temporary z positions
        Centroid_MID_cell_LUMP_TMP{c_f} = Centroid_MID_cell_LUMP_ORIG_UP{query_var_over_segTMP};
    end
    
% CASE 2: Only query fiber exists
elseif ~isempty(query_var) &&  isempty(query_var_over_segTMP)
    
    [NEW_CENTROID] = probe_var_correction_up(Centroid_MID_cell_LUMP_ORIG_UP,c_f,query_var,main_res,z_TMP);
    c_f = query_var;
    
% CASE 3: Only overlapping fiber exists
elseif isempty(query_var) &&  ~isempty(query_var_over_segTMP)
    
    % Define angle thresholds based on slice size
    Ang_query_slice_thresh = 10;
    ang_thresh_mag = 35;
    
    if size(Centroid_MID_cell_LUMP_TMP{query_var},1)<Ang_query_slice_thresh
        angle_thresh_tmp = ang_thresh_mag;
    else
        angle_thresh_tmp = ang_thresh;
    end    
    
    % Check if angular difference is within threshold
    if d_vec_angles_overTMP <= angle_thresh_tmp
        [NEW_CENTROID,Centroid_MID_cell_LUMP_ORIG_UP{query_var_over_segTMP},...
            z{query_var_over_segTMP},STATS_Pix_LOC_LUMP_UP{query_var_over_segTMP}] = Sticth_Shared_fibers(Centroid_MID_cell_LUMP_ORIG_UP,z,c_f,...
                                                                                                          query_var_over_segTMP,STATS_Pix_LOC_LUMP_UP);
        c_f = query_var_over_segTMP;
        z_TMP{c_f} = z{query_var_over_segTMP};      % Update temporary z positions
        Centroid_MID_cell_LUMP_TMP{c_f} = Centroid_MID_cell_LUMP_ORIG_UP{query_var_over_segTMP};
    else
        NEW_CENTROID = Centroid_MID_cell_LUMP_ORIG_UP{c_f};
        c_f = [];
    end
    
% CASE 4: Neither query nor overlapping fiber exists
elseif isempty(query_var) &&  isempty(query_var_over_segTMP)
    
    NEW_CENTROID = Centroid_MID_cell_LUMP_ORIG_UP{c_f};
    c_f = query_var;
end

% Return updated fiber information
varargout{1} = NEW_CENTROID;
varargout{end+1} = c_f;
varargout{end+1} = z;
varargout{end+1} = z_TMP;
varargout{end+1} = Centroid_MID_cell_LUMP_ORIG_UP;
varargout{end+1} = Centroid_MID_cell_LUMP_TMP;
varargout{end+1} = STATS_Pix_LOC_LUMP_UP;