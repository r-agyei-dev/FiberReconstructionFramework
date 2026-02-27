                                                                                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                                                            %% ALGORITHM FOR 3D RECONSTRUCTION CODE:
                                                                                                            % FINAL VERSION UPDATE: 11/12/2018.
                                                                                                            % BY: Ronald F Agyei
                                                                                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% This mfile performs the split of reconstructed fibers based on the overlap of the linear indices of the constituent ellipses 

function [varargout] = SPLITTER_CORRECTOR_FUNC_CODE(varargin)

Linear_Index_center_GLOBAL                 =      varargin{1};
minor_axis_length_extract                  =      varargin{2};
ORIGI_Centroid_MID_cel                     =      varargin{3};
height_variable                            =      varargin{4};
Pixel_Idx_Per_Slice_cell_GLOBAL            =      varargin{5};

% ================= WAITBAR INITIALIZATION =================
wb = waitbar(0,'Initializing Splitter Corrector...','Name','3D Reconstruction Progress');
set(wb,'Units','normalized','Position',[0.35 0.45 0.40 0.08]);
updateWB = @(frac,msg) waitbar(min(max(frac,0),1), wb, msg);
totalLoops = size(Linear_Index_center_GLOBAL,2);
% ===========================================================

Linear_Index_center_GLOBAL_NEW             =      cell(1,size(Linear_Index_center_GLOBAL,2));
ORIGI_Centroid_MID_cel_NEW                 =      cell(1,size(Linear_Index_center_GLOBAL,2));
Pixel_Idx_Per_Slice_cell_GLOBAL_NEW        =      cell(1,size(Linear_Index_center_GLOBAL,2));
numel_fib                                  =      cell(1,size(Linear_Index_center_GLOBAL,2));
minor_axis_length_extract_NEW              =      cell(1,size(Linear_Index_center_GLOBAL,2));

% loop through the variable of the linear indices of the fibers 
for ii = 1:size(Linear_Index_center_GLOBAL,2)
    
    % -------- WAITBAR UPDATE --------
    updateWB(ii/totalLoops, sprintf('Processing fiber %d of %d', ii, totalLoops));
    % --------------------------------

%    sprintf('%s%d','BEGIN Splitter corrector CUM defective fiber removal RUN ', ii)

        if ~iscell(Linear_Index_center_GLOBAL{ii})

            % Create temporary variables for the
            % (a) Linear Indices_variables         
            % (b) height_vec_variables 
            % (c) minor_vec_variables 
            % (d) ORIGI_Centroid_MID_cel_variables
            % (e) Unique_height_vec variables
                        
                                    if ~isempty(varargin{2}) 
                                    minor_vec_tmp                 =     minor_axis_length_extract{ii};
                                    end
            
                                    linear_vec_tmp                =     Linear_Index_center_GLOBAL{ii}; 
                                    height_vec_tmp                =     height_variable{ii};
                                    ORIGI_Centroid_MID_cel_tmp    =     ORIGI_Centroid_MID_cel{ii};
                                    unique_height_vec             =     unique(height_vec_tmp);            
                                    pptt                          =     Pixel_Idx_Per_Slice_cell_GLOBAL{ii};
                                    
                                    if ~iscell(pptt)
                                        pptt = {pptt};
                                    end
                                    
                                    if all(cell2mat(cellfun(@(x)isempty(x),pptt,'uni',0)) == 0)                   
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                
                                            intersect_var                 =     cell2mat(cellfun(@(x,y)~isempty(intersect(x,y)),pptt(1:end-1),pptt(2:end),'uni',0));

                                            % check to see if its an empty set or not

                                            numel_idx = 1:numel(Pixel_Idx_Per_Slice_cell_GLOBAL{ii});
                                            if any(intersect_var == 0)
                                                  [chunks_fibers_idx_cell_1] =  zeros_function_locater(intersect_var); 
                                                   chunks_fibers_idx_cell_1  =  cellfun(@(x)[x x(end)+1],chunks_fibers_idx_cell_1,'uni',0);  
                                                   chunks_fibers_idx_cell_2  =  [chunks_fibers_idx_cell_1    num2cell(numel_idx(~ismember(numel_idx,[chunks_fibers_idx_cell_1{:}])),1)];  
                                                   chunks_fibers_idx_cell_2  =  cellfun(@(x)unique_height_vec(x),chunks_fibers_idx_cell_2,'uni',0); 
                                                   
                                                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                   if ii == 1              
                                                   counter_vec = 1:numel(chunks_fibers_idx_cell_2);                   
                                                   else                    
                                                   counter_vec = counter_end + (1:numel(chunks_fibers_idx_cell_2));   
                                                   end
                                                   [numel_fib{ii}] = counter_vec;
                                                   counter_end = counter_vec(end);
                                                   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                      
                                                   Linear_Index_center_GLOBAL_NEW{ii}        =     cellfun(@(x)linear_vec_tmp(ismember(height_vec_tmp,x)),chunks_fibers_idx_cell_2,'uni',0);

                                                   if ~isempty(varargin{2})
                                                   minor_axis_length_extract_NEW{ii}         =     cellfun(@(x)minor_vec_tmp(ismember(unique_height_vec,x)),chunks_fibers_idx_cell_2,'uni',0);
                                                   end
                                                   
                                                   ORIGI_Centroid_MID_cel_NEW{ii}            =     cellfun(@(x)ORIGI_Centroid_MID_cel_tmp(ismember(unique_height_vec,x),:),chunks_fibers_idx_cell_2,'uni',0); 
                                                   
                                                   Pixel_Idx_Per_Slice_cell_GLOBAL_NEW{ii}   =     cellfun(@(x)pptt(ismember(unique_height_vec,x)),chunks_fibers_idx_cell_2,'uni',0);

                                            else 

                                                Linear_Index_center_GLOBAL_NEW{ii}           =     linear_vec_tmp;

                                                if ~isempty(varargin{2}) 
                                                minor_axis_length_extract_NEW{ii}            =     minor_vec_tmp;                                     
                                                end
                                                
                                                ORIGI_Centroid_MID_cel_NEW{ii}               =     ORIGI_Centroid_MID_cel_tmp;     
                                                Pixel_Idx_Per_Slice_cell_GLOBAL_NEW{ii}      =     pptt;
                                                
                                                if ii == 1
                                                counter_vec = 1;                
                                                else 
                                                counter_vec = counter_end + 1;  
                                                end
                                                [numel_fib{ii}] = counter_vec;
                                                counter_end = counter_vec(end);                                                
                                            end
                                            
                                    end
                      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
                                                                        
        else 
                    def_idx = ones(1,numel(Linear_Index_center_GLOBAL{ii}));

                    for jj = 1:numel(Linear_Index_center_GLOBAL{ii})   
                        if all(cell2mat(cellfun(@(x)isempty(x),Pixel_Idx_Per_Slice_cell_GLOBAL{ii}{jj},'uni',0)) == 0)      
                           def_idx(jj) = 0;      
                        end
                    end                    
                    
                    if all(def_idx == 0)

                                for jj = 1:numel(Linear_Index_center_GLOBAL{ii})

                                                linear_vec_tmp                =     Linear_Index_center_GLOBAL{ii}{jj};
                                                height_vec_tmp                =     height_variable{ii}{jj};
                                                
                                                if ~isempty(varargin{2}) 
                                                minor_vec_tmp                 =     minor_axis_length_extract{ii}{jj};
                                                end
                                                
                                                ORIGI_Centroid_MID_cel_tmp    =     ORIGI_Centroid_MID_cel{ii}{jj};
                                                unique_height_vec             =     unique(height_vec_tmp);
                                                pptt                          =     Pixel_Idx_Per_Slice_cell_GLOBAL{ii}{jj};

                                                intersect_var                 =     cell2mat(cellfun(@(x,y)~isempty(intersect(x,y)),pptt(1:end-1),pptt(2:end),'uni',0));

                                                numel_idx = 1:numel(Pixel_Idx_Per_Slice_cell_GLOBAL{ii}{jj});
                                                if any(intersect_var == 0)
                                                    [chunks_fibers_idx_cell_1] =  zeros_function_locater(intersect_var);
                                                    chunks_fibers_idx_cell_1  =  cellfun(@(x)[x x(end)+1],chunks_fibers_idx_cell_1,'uni',0);
                                                    chunks_fibers_idx_cell_2  =  [chunks_fibers_idx_cell_1    num2cell(numel_idx(~ismember(numel_idx,[chunks_fibers_idx_cell_1{:}])),1)];
                                                    chunks_fibers_idx_cell_2  =  cellfun(@(x)unique_height_vec(x),chunks_fibers_idx_cell_2,'uni',0);
                                                    
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
                                                    
                                                    Linear_Index_center_GLOBAL_NEW{ii}{jj}         =       cellfun(@(x)linear_vec_tmp(ismember(height_vec_tmp,x)),chunks_fibers_idx_cell_2,'uni',0);

                                                   if ~isempty(varargin{2})
                                                    minor_axis_length_extract_NEW{ii}{jj}          =       cellfun(@(x)minor_vec_tmp(ismember(unique_height_vec,x)),chunks_fibers_idx_cell_2,'uni',0);
                                                   end
                                                   
                                                    ORIGI_Centroid_MID_cel_NEW{ii}{jj}             =       cellfun(@(x)ORIGI_Centroid_MID_cel_tmp(ismember(unique_height_vec,x),:),chunks_fibers_idx_cell_2,'uni',0);
                                                    Pixel_Idx_Per_Slice_cell_GLOBAL_NEW{ii}{jj}     =       cellfun(@(x)pptt(ismember(unique_height_vec,x)),chunks_fibers_idx_cell_2,'uni',0);
                                                   
                                                else 
                                                    Linear_Index_center_GLOBAL_NEW{ii}{jj}    =     linear_vec_tmp;

                                                    if ~isempty(varargin{2})
                                                    minor_axis_length_extract_NEW{ii}{jj}     =     minor_vec_tmp;
                                                    end
                                                    
                                                    ORIGI_Centroid_MID_cel_NEW{ii}{jj}        =     ORIGI_Centroid_MID_cel_tmp;
                                                    Pixel_Idx_Per_Slice_cell_GLOBAL_NEW{ii}{jj}       =      pptt;  
                                                    
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
        
clearvars    chunks_fibers_idx_cell_1     chunks_fibers_idx_cell_2       linear_vec_tmp                     height_vec_tmp 
clearvars    minor_vec_tmp                Centroid_MID_cell_tmp          ORIGI_Centroid_MID_cel_tmp         unique_height_vec
clearvars    pptt                         intersect_var             

% sprintf('%s%d','END Splitter corrector CUM defective fiber removal RUN ', ii)

end

% -------- CLOSE WAITBAR --------
if isvalid(wb)
    close(wb);
end
% --------------------------------

varargout{            1}       =     Linear_Index_center_GLOBAL_NEW;
varargout{        end+1}       =     minor_axis_length_extract_NEW;
varargout{        end+1}       =     ORIGI_Centroid_MID_cel_NEW;
varargout{        end+1}       =     Pixel_Idx_Per_Slice_cell_GLOBAL_NEW;
varargout{        end+1}       =     numel_fib;

end