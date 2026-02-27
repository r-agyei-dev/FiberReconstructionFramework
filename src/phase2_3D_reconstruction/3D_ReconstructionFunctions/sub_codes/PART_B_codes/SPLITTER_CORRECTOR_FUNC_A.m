                                                                                                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                                                           %% SPLITTER_CORRECTOR_FUNC_A
                                                                                                           % 3D Fiber Reconstruction Splitter & Defective Fiber Corrector
                                                                                                           % FINAL VERSION UPDATE: 11/12/2018
                                                                                                           % BY: Ronald F Agyei
                                                                                                           %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This function splits reconstructed fibers based on overlaps in their
% linear indices and removes defective fibers. It ensures that fibers
% with non-overlapping segments are treated as separate fibers in the
% reconstruction pipeline.

function [varargout] = SPLITTER_CORRECTOR_FUNC_A(varargin)

% Input unpacking
Linear_Index_center_GLOBAL          = varargin{1};   % Original linear indices of fiber pixels
minor_axis_length_extract           = varargin{2};   % Minor axis lengths of ellipses
ORIGI_Centroid_MID_cel              = varargin{3};   % Original centroids of fiber ellipses
height_variable                     = varargin{4};   % Slice indices (height) of fiber pixels
Pixel_Idx_Per_Slice_cell_GLOBAL     = varargin{5};   % Pixel indices per slice for each fiber

% Initialize outputs with same structure as input
Linear_Index_center_GLOBAL_NEW       = cell(1,size(Linear_Index_center_GLOBAL,2));
minor_axis_length_extract_NEW        = cell(1,size(Linear_Index_center_GLOBAL,2));
ORIGI_Centroid_MID_cel_NEW           = cell(1,size(Linear_Index_center_GLOBAL,2));
Pixel_Idx_Per_Slice_cell_GLOBAL_NEW  = cell(1,size(Linear_Index_center_GLOBAL,2));
numel_fib                            = cell(1,size(Linear_Index_center_GLOBAL,2));

% Loop through each fiber's linear indices
for ii = 1:size(Linear_Index_center_GLOBAL,2)
    
    % Check if the fiber data is not a cell array (single fiber)
    if ~iscell(Linear_Index_center_GLOBAL{ii})

        % Assign temporary variables for this fiber
        linear_vec_tmp              = Linear_Index_center_GLOBAL{ii}; 
        height_vec_tmp              = height_variable{ii};
        minor_vec_tmp               = minor_axis_length_extract{ii};
        ORIGI_Centroid_MID_cel_tmp  = ORIGI_Centroid_MID_cel{ii};
        unique_height_vec           = unique(height_vec_tmp);            
        pptt                        = Pixel_Idx_Per_Slice_cell_GLOBAL{ii};
        
        % Only process if the fiber has non-empty pixel indices
        if all(cell2mat(cellfun(@(x)isempty(x),pptt,'uni',0)) == 0)                   
            % Determine overlapping segments between consecutive slices
            intersect_var = cell2mat(cellfun(@(x,y)~isempty(intersect(x,y)),pptt(1:end-1),pptt(2:end),'uni',0));
            
            numel_idx = 1:numel(Pixel_Idx_Per_Slice_cell_GLOBAL{ii});
            
            if any(intersect_var == 0)
                % Identify disconnected segments within the fiber
                [chunks_fibers_idx_cell_1] = zeros_function_locater(intersect_var); 
                chunks_fibers_idx_cell_1  = cellfun(@(x)[x x(end)+1],chunks_fibers_idx_cell_1,'uni',0);  
                chunks_fibers_idx_cell_2  = [chunks_fibers_idx_cell_1 num2cell(numel_idx(~ismember(numel_idx,[chunks_fibers_idx_cell_1{:}])),1)];  
                chunks_fibers_idx_cell_2  = cellfun(@(x)unique_height_vec(x),chunks_fibers_idx_cell_2,'uni',0); 
                
                % Determine counter vector for fiber numbering
                if ii == 1              
                    counter_vec = 1:numel(chunks_fibers_idx_cell_2);  
                else                    
                    counter_vec = counter_end + (1:numel(chunks_fibers_idx_cell_2));  
                end
                [numel_fib{ii}] = counter_vec;
                counter_end = counter_vec(end);
                
                % Split fiber variables into chunks based on overlaps
                Linear_Index_center_GLOBAL_NEW{ii}      = cellfun(@(x)linear_vec_tmp(ismember(height_vec_tmp,x)),chunks_fibers_idx_cell_2,'uni',0);
                minor_axis_length_extract_NEW{ii}       = cellfun(@(x)minor_vec_tmp(ismember(unique_height_vec,x)),chunks_fibers_idx_cell_2,'uni',0);
                ORIGI_Centroid_MID_cel_NEW{ii}         = cellfun(@(x)ORIGI_Centroid_MID_cel_tmp(ismember(unique_height_vec,x),:),chunks_fibers_idx_cell_2,'uni',0); 
                Pixel_Idx_Per_Slice_cell_GLOBAL_NEW{ii}= cellfun(@(x)pptt(ismember(unique_height_vec,x)),chunks_fibers_idx_cell_2,'uni',0);

            else
                % No overlaps detected, keep fiber as is
                Linear_Index_center_GLOBAL_NEW{ii}        = linear_vec_tmp;
                minor_axis_length_extract_NEW{ii}         = minor_vec_tmp;                                     
                ORIGI_Centroid_MID_cel_NEW{ii}           = ORIGI_Centroid_MID_cel_tmp;     
                Pixel_Idx_Per_Slice_cell_GLOBAL_NEW{ii}  = pptt;
                
                if ii == 1
                    counter_vec = 1;                 
                else 
                    counter_vec = counter_end + 1;   
                end
                [numel_fib{ii}] = counter_vec;
                counter_end = counter_vec(end); 
            end
        end
        
    else
        % Case: fiber data is already a cell array (multi-segment fiber)
        def_idx = ones(1,numel(Linear_Index_center_GLOBAL{ii}));

        % Identify defective sub-fibers
        for jj = 1:numel(Linear_Index_center_GLOBAL{ii})   
            if all(cell2mat(cellfun(@(x)isempty(x),Pixel_Idx_Per_Slice_cell_GLOBAL{ii}{jj},'uni',0)) == 0)      
               def_idx(jj) = 0;      
            end
        end                    

        % Only process if all sub-fibers are non-defective
        if all(def_idx == 0)
            for jj = 1:numel(Linear_Index_center_GLOBAL{ii})

                linear_vec_tmp               = Linear_Index_center_GLOBAL{ii}{jj};
                height_vec_tmp               = height_variable{ii}{jj};
                minor_vec_tmp                = minor_axis_length_extract{ii}{jj};
                ORIGI_Centroid_MID_cel_tmp   = ORIGI_Centroid_MID_cel{ii}{jj};
                unique_height_vec            = unique(height_vec_tmp);
                pptt                         = Pixel_Idx_Per_Slice_cell_GLOBAL{ii}{jj};

                intersect_var = cell2mat(cellfun(@(x,y)~isempty(intersect(x,y)),pptt(1:end-1),pptt(2:end),'uni',0));
                numel_idx = 1:numel(Pixel_Idx_Per_Slice_cell_GLOBAL{ii}{jj});

                if any(intersect_var == 0)
                    % Split the multi-segment fiber into disconnected chunks
                    [chunks_fibers_idx_cell_1] = zeros_function_locater(intersect_var);
                    chunks_fibers_idx_cell_1  = cellfun(@(x)[x x(end)+1],chunks_fibers_idx_cell_1,'uni',0);
                    chunks_fibers_idx_cell_2  = [chunks_fibers_idx_cell_1 num2cell(numel_idx(~ismember(numel_idx,[chunks_fibers_idx_cell_1{:}])),1)];
                    chunks_fibers_idx_cell_2  = cellfun(@(x)unique_height_vec(x),chunks_fibers_idx_cell_2,'uni',0);

                    % Fiber numbering counter
                    if ii == 1
                        if jj == 1
                            counter_vec = 1:numel(chunks_fibers_idx_cell_2);  
                        else
                            counter_vec = counter_end + (1:numel(chunks_fibers_idx_cell_2));  
                        end
                    else
                        counter_vec = counter_end + (1:numel(chunks_fibers_idx_cell_2));     
                    end
                    [numel_fib{ii}{jj}] = counter_vec;

                    % Assign split fiber data
                    Linear_Index_center_GLOBAL_NEW{ii}{jj}       = cellfun(@(x)linear_vec_tmp(ismember(height_vec_tmp,x)),chunks_fibers_idx_cell_2,'uni',0);
                    minor_axis_length_extract_NEW{ii}{jj}        = cellfun(@(x)minor_vec_tmp(ismember(unique_height_vec,x)),chunks_fibers_idx_cell_2,'uni',0);
                    ORIGI_Centroid_MID_cel_NEW{ii}{jj}          = cellfun(@(x)ORIGI_Centroid_MID_cel_tmp(ismember(unique_height_vec,x),:),chunks_fibers_idx_cell_2,'uni',0);
                    Pixel_Idx_Per_Slice_cell_GLOBAL_NEW{ii}{jj} = cellfun(@(x)pptt(ismember(unique_height_vec,x)),chunks_fibers_idx_cell_2,'uni',0);

                else
                    % No splits needed, keep fiber as is
                    Linear_Index_center_GLOBAL_NEW{ii}{jj}       = linear_vec_tmp;
                    minor_axis_length_extract_NEW{ii}{jj}        = minor_vec_tmp;
                    ORIGI_Centroid_MID_cel_NEW{ii}{jj}          = ORIGI_Centroid_MID_cel_tmp;
                    Pixel_Idx_Per_Slice_cell_GLOBAL_NEW{ii}{jj} = pptt;  

                    if ii == 1
                        if jj == 1
                            counter_vec = 1;
                        else
                            counter_vec = counter_end + 1; 
                        end
                    else 
                        counter_vec = counter_end + 1;    
                    end
                    [numel_fib{ii}{jj}] = counter_vec;
                end 
                counter_end = counter_vec(end);
                clearvars counter_vec
            end
        end
    end

    % Clear temporary variables to save memory
    clearvars chunks_fibers_idx_cell_1 chunks_fibers_idx_cell_2 linear_vec_tmp height_vec_tmp 
    clearvars minor_vec_tmp ORIGI_Centroid_MID_cel_tmp unique_height_vec pptt intersect_var

end

% Assign outputs
varargout{1} = Linear_Index_center_GLOBAL_NEW;
varargout{end+1} = minor_axis_length_extract_NEW;
varargout{end+1} = ORIGI_Centroid_MID_cel_NEW;
varargout{end+1} = Pixel_Idx_Per_Slice_cell_GLOBAL_NEW;
varargout{end+1} = numel_fib;
end