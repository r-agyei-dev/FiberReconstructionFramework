% ST_QUERY_FURTHER_ANALYSIS_F
% This function performs a refined analysis of a query fiber in a 3D fiber
% reconstruction context, similar to the '_b' variant. It focuses on determining
% the correct query fiber based on slice overlap, 3D orientation, in-plane
% distance, fiber length, and uniqueness criteria. This version considers
% the bottom-most slice of the fiber for overlap checks.
%
% Inputs:
%   query_var                       - Index of candidate query fiber(s)
%   crude_bin_mat                    - Flag indicating if overlap-based selection is used
%   Linear_Index_center_GLOBAL_LUMP_TMP - Linear indices of candidate fibers
%   size_length                      - Dimensions of the 3D image volume
%   main_res                         - Full set of main result coordinates
%   main_res_extract                 - Subset of main_res for the query fiber
%   rowsInImage, columnsInImage      - Image dimensions for 2D projections
%   radius                           - Radius used for in-plane distance checks
%   z_TMP                            - Z positions of fibers
%   Centroid_MID_cell_LUMP_TMP       - Centroid coordinates of candidate fibers
%   TEMP_MAT_NEW                     - Temporary matrix used in projection loops
%   Structure_Rogue                  - Information about rogue fibers
%   ORIGINAL_CRUDE_BIN_MAT           - Original binary projection matrix
%   direction_vectors                - 3D direction vectors for fibers
%   c_f                              - Current fiber index
%   inplane_thresh                   - Threshold for in-plane centroid distance
%   ang_thresh                        - Threshold for angular difference
%   circlePixels_points_IDX           - Linear indices for overlap check pixels

function [varargout] = ST_query_further_analysis_f(varargin)

query_var                                = varargin{1};
crude_bin_mat                            = varargin{2};
Linear_Index_center_GLOBAL_LUMP_TMP      = varargin{3};
size_length                              = varargin{4};
main_res                                 = varargin{5};
main_res_extract                         = varargin{6};
rowsInImage                              = varargin{7};
columnsInImage                           = varargin{8};
radius                                   = varargin{9};
z_TMP                                    = varargin{10};
Centroid_MID_cell_LUMP_TMP               = varargin{11};
TEMP_MAT_NEW                             = varargin{12};
Structure_Rogue                          = varargin{13};
ORIGINAL_CRUDE_BIN_MAT                   = varargin{14};
direction_vectors                        = varargin{15};
c_f                                      = varargin{16};
inplane_thresh                           = varargin{17};
ang_thresh                               = varargin{18}; 
circlePixels_points_IDX                  = varargin{19};

%% STEP 1: Update main_res and main_res_extract for query fiber
if ~isempty(query_var)
    
    if crude_bin_mat == 1
        % Determine the slice height based on overlap
        LIN_extract_lim = Linear_Index_center_GLOBAL_LUMP_TMP(query_var);
        [xt,yt,zt] = cellfun(@(x)ind2sub(size_length,x),LIN_extract_lim,'uni',0);
        clearvars xt yt
        
        % Use the first z-index of the fiber for overlap calculations
        query_height = cellfun(@(x)x(1), z_TMP(query_var), 'uni', 0);
        LIN_extract_lim = cellfun(@(x,y,z)y(ismember(x,z)), zt, LIN_extract_lim, query_height, 'uni',0);
        query_height = cell2mat(query_height);
        
        % Compute the overlap between projected fibers and main_res
        overlap_lim = zeros(1,numel(query_height));
        for qq = 1:numel(query_height)
            main_res_lim = main_res(main_res(:,3)== query_height(qq),1:2);
            if ~isempty(main_res_lim)
                Centers_pouplate_lim = find((columnsInImage - main_res_lim(1)).^2 + (rowsInImage - main_res_lim(2)).^2 <= radius.^2);
                clearvars circlePixels_cr_lim
                [circlePixels_cr_lim(:,2),circlePixels_cr_lim(:,1)] = ind2sub([size_length(1) size_length(2)],Centers_pouplate_lim);
                circlePixels_points_IDX_lim = sub2ind(size_length, circlePixels_cr_lim(:,2), circlePixels_cr_lim(:,1), query_height(qq)*(ones(size(circlePixels_cr_lim(:,1),1),1)));
                overlap_lim(qq) = numel(find(ismember(circlePixels_points_IDX_lim, LIN_extract_lim{qq})));
            else
                overlap_lim(qq) = -1;
            end
        end
        
        % Remove fibers whose projected slices do not overlap
        query_height(overlap_lim < 0) = [];
        query_var(overlap_lim < 0)    = [];
        overlap_lim(overlap_lim < 0)  = [];
        
        % Select the slice with maximum overlap and lowest height
        query_height_select = query_height(query_height == min(query_height(overlap_lim == max(overlap_lim))));
        query_var(query_height ~= query_height_select) = [];
        
        % Update main result matrices for the selected slice
        main_res_extract = main_res(main_res(:,3) == query_height_select,:);
        main_res(main_res(:,3) >= query_height_select,:) = [];
    end
end

%% STEP 2: Project through to seed fiber and extract coordinates
% Only perform projection if fiber is sufficiently long
Angle_query_slice_thresh = 10;
ang_thresh_mag = 35;

if ~isempty(query_var)
    zf = z_TMP{query_var};
    halt_var = 1;
    slice_limit = 0;
    query_path = 1;
    
    % Prepare container for projected main results
    clearvars main_res_extractV
    main_res_extractV = cell(1,numel(query_var));
    
    for gg = 1:numel(query_var)
        [~,main_res_extractN,main_resN,probe_varN,~, crude_bin_mat] = CENTERS_EXTRAPOLATE_LOOP_FUNC_NEW_UPDATE(halt_var,zf,slice_limit,query_var(gg),radius,columnsInImage, ...
            rowsInImage,size_length,Centroid_MID_cell_LUMP_TMP,query_path,TEMP_MAT_NEW,Structure_Rogue,z_TMP,ORIGINAL_CRUDE_BIN_MAT);
        
        if ~isempty(probe_varN) && crude_bin_mat == 1
            % Match z-position of the probe fiber
            main_res_extractN = main_resN(ismember(main_resN(:,3), z_TMP{probe_varN}(end)),:);
        else
            main_res_extractN = [];
        end
        
        main_res_extractV{gg} = main_res_extractN;
    end
    
    % Remove query fibers with empty projections
    query_var(cell2mat(cellfun(@(x)isempty(x), main_res_extractV, 'uni',0))) = [];
    main_res_extractV(cell2mat(cellfun(@(x)isempty(x), main_res_extractV, 'uni',0))) = [];
    
    %% STEP 2A: Angle measurement between seed and query fibers
    d_vec_cell = {direction_vectors{query_var}};
    d_vec_angles = cell2mat(cellfun(@(x)acosd(dot(direction_vectors{c_f},x)/(norm(direction_vectors{c_f})*norm(x))), d_vec_cell, 'uni',0));
    
    %% STEP 2B: In-plane centroid distance checks
    cent_cell_temp = {Centroid_MID_cell_LUMP_TMP{query_var}};
    dist_vec = cell2mat(cellfun(@(x)norm(main_res_extract(1,1:2) - x(1,1:2)), cent_cell_temp, 'uni',0));
    
    cent_cell_tempN = Centroid_MID_cell_LUMP_TMP{c_f};
    dist_vecN = cell2mat(cellfun(@(x)norm(x(1,1:2) - cent_cell_tempN(end,1:2)), main_res_extractV, 'uni',0));
    
    %% STEP 2C: Filter query fibers based on angular and in-plane thresholds
    if size(Centroid_MID_cell_LUMP_TMP{c_f},1) > Angle_query_slice_thresh
        query_var(~(d_vec_angles<ang_thresh & (dist_vec<inplane_thresh | dist_vecN<inplane_thresh))) = [];
    else
        query_var(~(d_vec_angles<ang_thresh_mag & (dist_vec<inplane_thresh | dist_vecN<inplane_thresh))) = [];
    end
end

%% STEP 3: Overlap check with seed fiber
if ~isempty(query_var)
    LIN_extract_tmp = Linear_Index_center_GLOBAL_LUMP_TMP(query_var);
    [xs,ys,zs] = cellfun(@(x)ind2sub(size_length,x), LIN_extract_tmp, 'uni',0);
    clearvars xs ys
    LIN_extract = cellfun(@(x,y)y(ismember(x,main_res_extract(end))), zs, LIN_extract_tmp, 'uni',0);
    LIN_extract_overlap = cell2mat(cellfun(@(x)numel(find(ismember(x,circlePixels_points_IDX)))/numel(circlePixels_points_IDX), LIN_extract,'uni',0));
    query_var(LIN_extract_overlap ~= max(LIN_extract_overlap)) = [];
end

%% STEP 4: Ensure uniqueness of the query fiber
% If multiple candidates, choose based on angular difference
if numel(query_var) > 1
    d_vec_cell = {direction_vectors{query_var}};
    d_vec_angles = cell2mat(cellfun(@(x)acosd(dot(direction_vectors{c_f},x)/(norm(direction_vectors{c_f})*norm(x))), d_vec_cell, 'uni',0));
    query_var(d_vec_angles == max(d_vec_angles)) = [];
end

% If still multiple, choose tallest fiber
if numel(query_var) > 1
    d_vec_cell = {Centroid_MID_cell_LUMP_TMP{query_var}};
    d_vec_lsizes = cell2mat(cellfun(@(x)size(x,1), d_vec_cell, 'uni',0));
    query_var(d_vec_lsizes == max(d_vec_lsizes)) = [];
end

% If still multiple, pick the first one
if numel(query_var) > 1
    query_var(2:end) = [];
end

%% STEP 5: Final Variable Output
varargout{1} = query_var;
varargout{end+1} = main_res_extract;
varargout{end+1} = main_res;