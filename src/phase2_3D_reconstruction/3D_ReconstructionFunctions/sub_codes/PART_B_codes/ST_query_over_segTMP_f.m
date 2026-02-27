% ST_QUERY_OVER_SEGTMP_F
% This function identifies candidate query fibers that may be affected by
% over-segmentation or encasing relative to a seed fiber (c_f) in the "forward" direction.
% It handles cases where fibers fully encapsulate each other, share slices,
% or partially overlap, and filters candidates based on angular alignment.
%
% Inputs:
%   c_f                               - Index of the current (seed) fiber
%   query_var                         - Candidate query fiber indices
%   Crude_Idx_Split_Idx_var            - Mapping of fibers to unsplit microstructure indices
%   z_TMP                              - Cell array of z-coordinates for each fiber
%   ze                                 - Linear indices corresponding to z-slices
%   Linear_Index_center_GLOBAL_LUMP_TMP - Linear indices of all fiber voxels
%   STATS_Pix_LOC_LUMP_UP              - Pixel-level stats for each fiber
%   direction_vectors                  - 3D orientation vectors for fibers
%   size_length                        - Size of the 3D volume

function [varargout] = ST_query_over_segTMP_f(varargin)

c_f                                   = varargin{1};
query_var                             = varargin{2};
Crude_Idx_Split_Idx_var               = varargin{3};
z_TMP                                 = varargin{4};
ze                                    = varargin{5};
Linear_Index_center_GLOBAL_LUMP_TMP   = varargin{6};
STATS_Pix_LOC_LUMP_UP                 = varargin{7};
direction_vectors                     = varargin{8};
size_length                           = varargin{9};

% Initialize vectors to track over-segmented and fully-encasing fibers
query_var_over_seg = zeros(numel(query_var),1);
query_encase_del   = zeros(numel(query_var),1);

for pp = 1:numel(query_var)
    
    query_var_TMP = query_var(pp);
    
    % Map probe and seed fibers to their original unsplit structure
    Crude_idx_probe = cell2mat(Crude_Idx_Split_Idx_var(1,cell2mat(cellfun(@(x)any(ismember(x,query_var_TMP)), Crude_Idx_Split_Idx_var(end,:), 'uni',0))));
    Crude_idx_c_f   = cell2mat(Crude_Idx_Split_Idx_var(1,cell2mat(cellfun(@(x)any(ismember(x,c_f)), Crude_Idx_Split_Idx_var(end,:), 'uni',0))));
    
    % Handle fibers that fully encapsulate each other:
    % (a) query encases probe
    % (b) probe encases query
    if (z_TMP{c_f}(1) >= z_TMP{query_var_TMP}(1) && z_TMP{c_f}(end) <= z_TMP{query_var_TMP}(end)) || ...
       (z_TMP{c_f}(1) <= z_TMP{query_var_TMP}(1) && z_TMP{c_f}(end) >= z_TMP{query_var_TMP}(end))
        
        query_encase_del(pp) = query_var_TMP;
        
    else
        % Partial overlap case: determine if fibers share slices and classify
        overlap_depth = 7;
        if (z_TMP{query_var_TMP}(1) <= z_TMP{c_f}(end))
            
            if (abs(z_TMP{c_f}(end) - z_TMP{query_var_TMP}(1)) + 1 <= overlap_depth)
                
                if Crude_idx_probe == Crude_idx_c_f
                    % Same unsplit fiber: mark for deletion
                    query_encase_del(pp) = query_var_TMP;
                else
                    overseg_thresh = 0.001;
                    
                    % Extract voxel indices for the first slice of the query fiber
                    LIN_extract_tmp_q = Linear_Index_center_GLOBAL_LUMP_TMP{query_var_TMP};
                    LIN_extract_tmp_q(ze{query_var_TMP} ~= z_TMP{query_var_TMP}(1)) = [];
                    
                    % Project probe fiber to query fiber start slice and find overlap
                    slice_loc = z_TMP{query_var_TMP}(1);
                    if nonzeros(ismember(z_TMP{c_f}(1):z_TMP{c_f}(end), slice_loc-1))
                        probe_extract_pix = STATS_Pix_LOC_LUMP_UP{c_f}(ismember(z_TMP{c_f}(1):z_TMP{c_f}(end), slice_loc-1)).FINAL_VALUES;
                        LIN_extract_tmp_p_proj = sub2ind(size_length, probe_extract_pix(:,2), probe_extract_pix(:,1), slice_loc*ones(size(probe_extract_pix,1),1));
                        
                        % Classify fiber based on projected overlap ratio
                        if numel(find(ismember(LIN_extract_tmp_q, LIN_extract_tmp_p_proj)))/numel(LIN_extract_tmp_p_proj) <= overseg_thresh
                            query_encase_del(pp) = query_var_TMP;
                        else
                            query_var_over_seg(pp) = query_var_TMP;
                        end
                    end
                end
                
            else
                % Excess overlap: mark for deletion
                query_encase_del(pp) = query_var_TMP;
            end
        else
            query_encase_del(pp) = 0;
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove fibers flagged for encasing or over-segmentation from the candidate set
query_var(ismember(query_var, query_encase_del))    = [];
query_var(ismember(query_var, query_var_over_seg))  = [];

% Select the fiber with the smallest out-of-plane angular difference
query_var_over_seg = nonzeros(query_var_over_seg);
if ~isempty(query_var_over_seg)
    d_vec_cell_over       = {direction_vectors{query_var_over_seg}};
    d_vec_angles_over     = cell2mat(cellfun(@(x)acosd(dot(direction_vectors{c_f},x)/(norm(direction_vectors{c_f})*norm(x))), d_vec_cell_over,'uni',0));
    query_var_over_segTMP = query_var_over_seg(d_vec_angles_over == min(d_vec_angles_over));
    d_vec_angles_overTMP  = d_vec_angles_over(d_vec_angles_over == min(d_vec_angles_over));
    
    % If multiple candidates remain, use overlap with TEMP_MAT_NEW to select best
    if numel(query_var_over_segTMP) > 1
        overlap = zeros(1,numel(query_var_over_segTMP));
        for nn = 1:numel(query_var_over_segTMP)
            overlap(nn) = numel(find(ismember(nonzeros(TEMP_MAT_NEW(circlePixels_points_IDX)), query_var_over_segTMP(nn))));
        end
        query_var_over_segTMP(overlap ~= max(overlap)) = [];
        d_vec_angles_overTMP(overlap ~= max(overlap))  = [];
    end
else
    query_var_over_segTMP = {};
    d_vec_angles_overTMP  = {};
end

%% Output variables
varargout{1} = query_var;
varargout{end+1} = query_var_over_segTMP;
varargout{end+1} = d_vec_angles_overTMP;