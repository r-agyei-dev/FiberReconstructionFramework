%% SPLITTER_CORRECTOR_FUNC_B
% This function corrects split fibers for reconstruction purposes. 
% It examines fibers represented by their linear indices, centroid positions, 
% and pixel-level information to detect discontinuities, then splits fibers
% into contiguous segments when necessary.
%
% Inputs:
%   Linear_Index_center_GLOBAL   - Cell array containing linear indices of fiber voxels
%   ORIGI_Centroid_MID_cel       - Cell array of original fiber centroid positions
%   Pixel_Idx_cell_GLOBAL        - Cell array of pixel indices for each fiber slice
%   height_variable              - Cell array of z-coordinates (slice positions) for fiber voxels
%
% Outputs:
%   varargout{1} - Corrected linear indices of fibers after splitting
%   varargout{2} - Corrected original centroids corresponding to split fibers

function [varargout] = SPLITTER_CORRECTOR_FUNC_B(varargin)

Linear_Index_center_GLOBAL       = varargin{1};
ORIGI_Centroid_MID_cel           = varargin{2};
Pixel_Idx_cell_GLOBAL            = varargin{3};
height_variable                  = varargin{4};

% Initialize output cells
Linear_Index_center_GLOBAL_NEW   = cell(1, size(Linear_Index_center_GLOBAL,2));
ORIGI_Centroid_MID_cel_NEW      = cell(1, size(Linear_Index_center_GLOBAL,2));

%% Iterate through all fibers
for ii = 1:size(Linear_Index_center_GLOBAL,2)
    
    % If fiber is not already a cell of sub-fibers
    if ~iscell(Linear_Index_center_GLOBAL{ii})
        
        % Extract relevant data for the current fiber
        linear_vec_tmp              = Linear_Index_center_GLOBAL{ii}; 
        height_vec_tmp              = height_variable{ii};
        ORIGI_Centroid_MID_cel_tmp = ORIGI_Centroid_MID_cel{ii};
        unique_height_vec           = unique(height_vec_tmp);            
        pptt                        = Pixel_Idx_cell_GLOBAL{ii};
        
        % Ensure pixel indices are in cell format
        if ~iscell(pptt)
            pptt = {pptt};
        end
        
        % Determine if there are non-overlapping slices between consecutive pixels
        intersect_var = cell2mat(cellfun(@(x,y)~isempty(intersect(x,y)), pptt(1:end-1), pptt(2:end), 'uni', 0));
        numel_idx = 1:numel(Pixel_Idx_cell_GLOBAL{ii});
        
        % If discontinuities exist, split the fiber into contiguous segments
        if any(intersect_var == 0)
            
            [chunks_fibers_idx_cell_1] = zeros_function_locater(intersect_var); 
            chunks_fibers_idx_cell_1   = cellfun(@(x)[x x(end)+1], chunks_fibers_idx_cell_1, 'uni', 0);  
            chunks_fibers_idx_cell_2   = [chunks_fibers_idx_cell_1 num2cell(numel_idx(~ismember(numel_idx,[chunks_fibers_idx_cell_1{:}])),1)];  
            chunks_fibers_idx_cell_2   = cellfun(@(x)unique_height_vec(x), chunks_fibers_idx_cell_2, 'uni', 0); 

            % Assign split segments to the outputs
            Linear_Index_center_GLOBAL_NEW{ii}   = cellfun(@(x)linear_vec_tmp(ismember(height_vec_tmp,x)), chunks_fibers_idx_cell_2, 'uni',0);
            ORIGI_Centroid_MID_cel_NEW{ii}      = cellfun(@(x)ORIGI_Centroid_MID_cel_tmp(ismember(unique_height_vec,x),:), chunks_fibers_idx_cell_2, 'uni',0);
        
        else
            % If no discontinuities, keep fiber as is
            Linear_Index_center_GLOBAL_NEW{ii}  = linear_vec_tmp;
            ORIGI_Centroid_MID_cel_NEW{ii}     = ORIGI_Centroid_MID_cel_tmp;     
        end
        
    else
        % If fiber is already subdivided into sub-cells
        for jj = 1:numel(Linear_Index_center_GLOBAL{ii})
            
            % Extract data for current sub-fiber
            linear_vec_tmp              = Linear_Index_center_GLOBAL{ii}{jj};
            height_vec_tmp              = height_variable{ii}{jj};
            ORIGI_Centroid_MID_cel_tmp = ORIGI_Centroid_MID_cel{ii}{jj};
            unique_height_vec           = unique(height_vec_tmp);
            pptt                        = Pixel_Idx_cell_GLOBAL{ii}{jj};
            
            if ~iscell(pptt)
                pptt = {pptt};
            end                                    
            
            % Check for slice discontinuities
            intersect_var = cell2mat(cellfun(@(x,y)~isempty(intersect(x,y)), pptt(1:end-1), pptt(2:end), 'uni',0));
            numel_idx = 1:numel(Pixel_Idx_cell_GLOBAL{ii}{jj});
            
            if any(intersect_var == 0)
                [chunks_fibers_idx_cell_1] = zeros_function_locater(intersect_var);
                chunks_fibers_idx_cell_1   = cellfun(@(x)[x x(end)+1], chunks_fibers_idx_cell_1,'uni',0);
                chunks_fibers_idx_cell_2   = [chunks_fibers_idx_cell_1 num2cell(numel_idx(~ismember(numel_idx,[chunks_fibers_idx_cell_1{:}])),1)];
                chunks_fibers_idx_cell_2   = cellfun(@(x)unique_height_vec(x), chunks_fibers_idx_cell_2,'uni',0);

                Linear_Index_center_GLOBAL_NEW{ii}{jj}   = cellfun(@(x)linear_vec_tmp(ismember(height_vec_tmp,x)), chunks_fibers_idx_cell_2,'uni',0);
                ORIGI_Centroid_MID_cel_NEW{ii}{jj}      = cellfun(@(x)ORIGI_Centroid_MID_cel_tmp(ismember(unique_height_vec,x),:), chunks_fibers_idx_cell_2,'uni',0);
            
            else
                % No discontinuity; keep sub-fiber as is
                Linear_Index_center_GLOBAL_NEW{ii}{jj}  = linear_vec_tmp;
                ORIGI_Centroid_MID_cel_NEW{ii}{jj}     = ORIGI_Centroid_MID_cel_tmp;
            end
        end
    end
    
    % Clear temporary variables for next iteration
    clearvars chunks_fibers_idx_cell_1 chunks_fibers_idx_cell_2 linear_vec_tmp height_vec_tmp 
    clearvars ORIGI_Centroid_MID_cel_tmp unique_height_vec pptt intersect_var
end

%% Assign outputs
varargout{1} = Linear_Index_center_GLOBAL_NEW;
varargout{2} = ORIGI_Centroid_MID_cel_NEW;

end