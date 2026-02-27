%% Version update: Changed the location where the refitted centroid is used (STATS_BREAK_FUNCTION_22_update)

% This function processes separated fiber stacks
% STEP 1: Further splits fibers based on out-of-plane orientation
% STEP 2: Corrects staggered centroids in 3D using SVD and Cook's distance
% STEP 3: Splits the structure file accordingly
% STEP 4: Fits ellipses around centroids to reconstruct 3D fiber volumes
% STEP 5: Generates a linear indexing for all fibers

%% CENTRAL_PIXEL_POPULATOR_FUNCTION
%   - Main orchestrator function for splitting fibers, correcting centroids,
%     extracting metadata, and populating pixels into 3D fiber structures.

function [varargout] = CENTRAL_PIXEL_POPULATOR_FUNCTION(varargin)

Distinct_Fiber_IDX_CELL_temp           = varargin{1};
GLOBAL_CENT_Distinct_Fiber_temp        = varargin{2};
main_cell                              = varargin{3};
fib_idx_n_slice_loc_cell               = varargin{4};
idx_after_first_slice_temp             = varargin{5};
size_length                            = varargin{6};
temp_thresh                            = varargin{7};

% disp('USING CENTRAL_PIXEL_GLOBAL_ESTIMATOR_UPDATED_DIRECT_DENSE_v2_TT.m FUNCTION')

% Copy initial data for further processing
fin_lump_VAL_CELL_TEMP                 = Distinct_Fiber_IDX_CELL_temp;
GLOBAL_CENT_TOTAL_TEMP                 = GLOBAL_CENT_Distinct_Fiber_temp;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Split fibers further based on out-of-plane orientation
% For each fiber, check if it can be subdivided using the 3D out-of-plane splitter.
% Outputs updated fin_lump_VAL_CELL_TEMP and GLOBAL_CENT_TOTAL_TEMP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clearvars bb Centroid_MID_cell_lump Pixel_units spacing spacing_depth moving_idx grad
diary('new_output.txt')

for index_cell = 1:size(Distinct_Fiber_IDX_CELL_temp,2)
    bb = GLOBAL_CENT_Distinct_Fiber_temp{index_cell};

    if iscell(bb)
        for gg = 1:length(bb)
            fig_switch = 0;
            [fin_lump_VAL_CELL_TEMP{index_cell}{gg}, GLOBAL_CENT_TOTAL_TEMP{index_cell}{gg}] = ...
                THREE_3D_OUT_PLANE_SPLITTER_DIRECT(GLOBAL_CENT_TOTAL_TEMP{index_cell}{gg}, fin_lump_VAL_CELL_TEMP{index_cell}{gg}, fig_switch);
        end
    else
        fig_switch = 0;
        [fin_lump_VAL_CELL_TEMP{index_cell}, GLOBAL_CENT_TOTAL_TEMP{index_cell}] = ...
            THREE_3D_OUT_PLANE_SPLITTER_DIRECT(GLOBAL_CENT_TOTAL_TEMP{index_cell}, fin_lump_VAL_CELL_TEMP{index_cell}, fig_switch); 
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Correct staggered centroids using 3D line fits (SVD + Cook's distance)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

GLOBAL_CENT_TOTAL_TEMP_refit = cell(1,length(GLOBAL_CENT_TOTAL_TEMP));
grad_out                     = cell(1,length(GLOBAL_CENT_TOTAL_TEMP));

for ii = 1:length(GLOBAL_CENT_TOTAL_TEMP)
    if iscell(GLOBAL_CENT_TOTAL_TEMP{ii})  
        for jj = 1:length(GLOBAL_CENT_TOTAL_TEMP{ii})
            if ~iscell(GLOBAL_CENT_TOTAL_TEMP{ii}{jj})
                [GLOBAL_CENT_TOTAL_TEMP_refit{ii}{jj}, grad_out{ii}{jj}] = ...
                    CENTROID_STAGGER_CORRRECTOR_DIRECT(GLOBAL_CENT_TOTAL_TEMP{ii}{jj}, fig_switch, size_length);
            else
                for hh = 1:length(GLOBAL_CENT_TOTAL_TEMP{ii}{jj})
                    [GLOBAL_CENT_TOTAL_TEMP_refit{ii}{jj}{hh}, grad_out{ii}{jj}{hh}] = ...
                        CENTROID_STAGGER_CORRRECTOR_DIRECT(GLOBAL_CENT_TOTAL_TEMP{ii}{jj}{hh}, fig_switch, size_length);
                end
            end
        end
    else 
        [GLOBAL_CENT_TOTAL_TEMP_refit{ii}, grad_out{ii}] = ...
            CENTROID_STAGGER_CORRRECTOR_DIRECT(GLOBAL_CENT_TOTAL_TEMP{ii}, fig_switch, size_length);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Prepare metadata structures and split the structure file
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

STATS_NEW_MIDDLE_LINE_GLOBAL    = cell(1,size(fin_lump_VAL_CELL_TEMP,2));
Linear_Index_center_GLOBAL      = cell(1,size(fin_lump_VAL_CELL_TEMP,2));
Pixel_Idx_Per_Slice_cell_LUMP   = cell(1,size(fin_lump_VAL_CELL_TEMP,2));
minRadius_idxval_LUMP           = cell(1,size(fin_lump_VAL_CELL_TEMP,2));
height_variable_cell            = cell(1,size(fin_lump_VAL_CELL_TEMP,2));

STATS_TEMP = cell(1,length(temp_thresh));
for ii = 1:size(idx_after_first_slice_temp,2)
    STATS_TEMP{ii} = main_cell{idx_after_first_slice_temp(1,ii)};
end

% Update fiber index mapping using STATS_BREAK_FUNCTION_22_update
[In_plane_angle_CELL_FINAL, fin_lump_VAL_CELL_TEMP_O, GLOBAL_CENT_TOTAL_TEMP_refit_UP, GLOBAL_CENT_TOTAL_TEMP_UP, ...
 grad_out_UP, Maj_Min_CELL_FINAL] = ...
    STATS_BREAK_FUNCTION_22_update(fin_lump_VAL_CELL_TEMP, STATS_TEMP, GLOBAL_CENT_TOTAL_TEMP_refit, GLOBAL_CENT_TOTAL_TEMP, grad_out);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: Fit ellipses around centroids and generate linear indexing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for mm = 1:size(fin_lump_VAL_CELL_TEMP_O,2)
    [iter_slice, ~] = RETRACKER_FUNCTION(fib_idx_n_slice_loc_cell, idx_after_first_slice_temp(mm));

    if iscell(fin_lump_VAL_CELL_TEMP_O{mm})
        STATS_NEW_MIDDLE_LINE = cell(1,length(fin_lump_VAL_CELL_TEMP_O{mm}));
        Pixel_Idx_cell        = cell(1,length(fin_lump_VAL_CELL_TEMP_O{mm}));
        minRadius_idxval      = cell(1,length(fin_lump_VAL_CELL_TEMP_O{mm}));
        defective             = zeros(1,length(fin_lump_VAL_CELL_TEMP_O{mm}));

        for ii = 1:length(fin_lump_VAL_CELL_TEMP_O{mm})
            GLOBAL_CENT_TOTAL_TEMP_UP{mm}{ii}(:,3) = GLOBAL_CENT_TOTAL_TEMP_refit_UP{mm}{ii}(:,3);
            [Pixel_Idx_cell{ii}, STATS_NEW_MIDDLE_LINE{ii}, minRadius_idxval{ii}, defective(ii)] = ...
                PIXEL_INTERPOLATE_FORMULATOR_SMART_ELLIPSE(In_plane_angle_CELL_FINAL{mm}{ii}, GLOBAL_CENT_TOTAL_TEMP_UP{mm}{ii}, size_length, Maj_Min_CELL_FINAL{mm}{ii});
        end
    else
        GLOBAL_CENT_TOTAL_TEMP_UP{mm}(:,3) = GLOBAL_CENT_TOTAL_TEMP_refit_UP{mm}(:,3);
        [Pixel_Idx_cell, STATS_NEW_MIDDLE_LINE, minRadius_idxval, defective] = ...
            PIXEL_INTERPOLATE_FORMULATOR_SMART_ELLIPSE(In_plane_angle_CELL_FINAL{mm}, GLOBAL_CENT_TOTAL_TEMP_UP{mm}, size_length, Maj_Min_CELL_FINAL{mm});
    end

    if all(defective==0)
        STATS_NEW_MIDDLE_LINE_GLOBAL{mm} = STATS_NEW_MIDDLE_LINE;
        Pixel_Idx_Per_Slice_cell_LUMP{mm} = Pixel_Idx_cell;
        minRadius_idxval_LUMP{mm} = minRadius_idxval;

        [check_mat_center, height_variable] = SPATIAL_INSTANCE_ASSEMBLER_AUTOMATOR_MULTI(STATS_NEW_MIDDLE_LINE, fin_lump_VAL_CELL_TEMP_O{mm}, iter_slice, size_length);

        Linear_Index_center_GLOBAL{mm} = check_mat_center;
        height_variable_cell{mm} = height_variable;
    end
end

% Output assignments
varargout{1} = STATS_NEW_MIDDLE_LINE_GLOBAL;
varargout{end+1} = Linear_Index_center_GLOBAL;
varargout{end+1} = GLOBAL_CENT_TOTAL_TEMP_refit_UP;
varargout{end+1} = GLOBAL_CENT_TOTAL_TEMP_UP;
varargout{end+1} = Pixel_Idx_Per_Slice_cell_LUMP;
varargout{end+1} = grad_out_UP;
varargout{end+1} = fin_lump_VAL_CELL_TEMP_O;
varargout{end+1} = minRadius_idxval_LUMP;
varargout{end+1} = height_variable_cell;


end
%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% PIPELINE GLOSSARY
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% GLOSSARY OF FUNCTIONS 
% Function 1:   ==>>   THREE_3D_OUT_PLANE_SPLITTER_DIRECT
% Function 2:   ==>>   CENTROID_STAGGER_CORRRECTOR_DIRECT
% Function 3:   ==>>   LINEAR_FIT_FUNCTION_3D_WITH_FIG_DIRECT
% Function 4:   ==>>   STATS_BREAK_FUNCTION_22_update
% Function 5:   ==>>   PIXEL_INTERPOLATE_FORMULATOR_SMART_ELLIPSE
% Function 6:   ==>>   SPATIAL_INSTANCE_ASSEMBLER_AUTOMATOR_MULTI
% Function 7:   ==>>   regstats_func
% Function 8:   ==>>   LIROWS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% THREE_3D_OUT_PLANE_SPLITTER_DIRECT
%   - Splits a fiber strand into sub-segments based on 3D orientation changes.
%   - Performs moving centroid linear fits to detect tortuous points.
%   - Outputs updated fiber indices and centroid coordinates.

function [varargout] = THREE_3D_OUT_PLANE_SPLITTER_DIRECT(varargin)

% This function evaluates the turtuosity of the strands of the fibers and perfroms the split.
% Using the aforementioned information it splits the fiber indices and the fiber centroids.

% Split the fibers into many parts and progressively fit lines through them
% and if any doesnt fit break

% The input for the function is extract of a strand of fiber centroids and further performs a seperation filter for the strand 
% Output can be a cell iplying the input strand was actually a a bifurcated fiber

GLOBAL_CENT_TOTAL_loop         =      varargin{1};
FIN_CENT_TOTAL_loop            =      varargin{2};
fig_switch                     =      varargin{3};

       % ====>>>>   STEP 1  : Declare the moving_idx variable 
      clearvars XY_NEW     moving_idx     grad_check         
      spacing                   =       10; 
      spacing_depth             =       spacing - 2;
      moving_idx                =       vertcat(1 : spacing_depth : size(GLOBAL_CENT_TOTAL_loop,1) - spacing , spacing + 1 : spacing_depth : size(GLOBAL_CENT_TOTAL_loop,1))';
            
  if ~isempty(moving_idx)   
      
      moving_idx(end)          =       size(GLOBAL_CENT_TOTAL_loop,1);        
      % ====>>>>   STEP 2  : Extract direction vectors for each of the segments                     
        direct_vect = zeros(size(moving_idx,1),3);
        for pp                 =          1        :       size(moving_idx,1)
        Pixel_units                                                                    =       GLOBAL_CENT_TOTAL_loop(moving_idx(pp,1):moving_idx(pp,end),:);
        [r_3D,~]                                                                       =       LINEAR_FIT_FUNCTION_3D_WITH_FIG_DIRECT(Pixel_units,fig_switch);
        direct_vect(pp,:)                                                              =       r_3D(1,:) - r_3D(end,:) ./ vecnorm(r_3D(1,:) - r_3D(end,:),2,2);  % normr(r_3D(1,:) - r_3D(end,:));
        end
        
     % ====>>>>   STEP 3    : Iteratively add up up the stacks using the direction vectors         
      clearvars fiber_break_keep 
      iter_val           =   1 : size(direct_vect,1)-1; % list all the numbers of the 
      break_count        =   0;
      keep               =   0;
      fiber_break_keep   =   cell(1,size(direct_vect,1));

       if isempty(iter_val)                              % for indexes where the moving index is of size 1
       FIN_CENT_TOTAL_loop_update      =   FIN_CENT_TOTAL_loop;
       updated_cent                    =   GLOBAL_CENT_TOTAL_loop;

       else
             while ~isempty(iter_val)
             pp               =   iter_val(1);
             count            =   1;
             keep(count)      =   pp;
             min_check        =   0;

                  while min_check < 45 && pp < size(direct_vect,1)
                  min_check       =     min([180 - acosd(dot(direct_vect(pp,:),direct_vect(pp+1,:)))  ;  acosd(dot(direct_vect(pp,:),direct_vect(pp+1,:))) ]);  
                  keep(count)     =     pp;
                  pp              =     pp + 1;
                  count           =     count +  1;
                  end
 
                  break_count                      =   break_count + 1;
                  fiber_break_keep{break_count}    =   keep;
                  iter_val(1:numel(keep))          =   [];
                  clearvars keep   count
            end

%  ====>>>>   STEP 5:  This part checks the end stock 
                  if dot(direct_vect(end,:),[0 0 1])== 0 ||  min_check < 30 
                     fiber_break_keep{break_count}(end + 1) = pp;
                  else 
                     fiber_break_keep{break_count+1} = pp;
                  end
                                    
% ====>>>>   STEP 6:   Delete the empty cells in the variable 
    fiber_break_keep(cell2mat(cellfun(@(x)isempty(x),fiber_break_keep,'uni',0))) = [];

% ====>>>>   STEP 7:   Extract the indices of the fibers after the splitting 
    split_fiber_idx = cellfun(@(x)x(1):x(end),cellfun(@(x)moving_idx(x,:),fiber_break_keep,'uni',0),'uni',0);

% ====>>>>   STEP 8:   Extract the indices of the fibers after the splitting 
    FIN_CENT_TOTAL_loop_update = cellfun(@(x)FIN_CENT_TOTAL_loop(x),split_fiber_idx,'uni',0);
    
% Update the centrids and the indices for the centroids for the fiber stalks
    updated_cent = cellfun(@(x)GLOBAL_CENT_TOTAL_loop(x,:),split_fiber_idx,'uni',0);

            if numel(FIN_CENT_TOTAL_loop_update)== 1
            FIN_CENT_TOTAL_loop_update  =  cell2mat(FIN_CENT_TOTAL_loop_update);
            updated_cent                =  cell2mat(updated_cent);
            end
    
       end
      
  else

    FIN_CENT_TOTAL_loop_update      =   FIN_CENT_TOTAL_loop;
    updated_cent                    =   GLOBAL_CENT_TOTAL_loop;

  end

varargout{               1}       =   FIN_CENT_TOTAL_loop_update;
varargout{           end+1}       =   updated_cent;


end



%% CENTROID_STAGGER_CORRRECTOR_DIRECT
%   - Corrects staggered centroids along fiber strands using SVD-based 3D line fits.
%   - Accounts for edge constraints and generates gradient information.

function [varargout] = CENTROID_STAGGER_CORRRECTOR_DIRECT(varargin)
    
    Centroid_MID_cell_lump   =   varargin{1};
    fig_switch               =   varargin{2};
    size_length              =   varargin{3};
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Please note that the ellipses need to be more than one otherwise the following error will be given by MATLAB
    % Error using regstats (line 124)
    % The design matrix has more predictor variables than observations.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if size(Centroid_MID_cell_lump,1) > 2
    % Perform the 3D fit to generate the points 
    
          Pixel_units_fin             =       Centroid_MID_cell_lump;
          [r_3D,grad]                 =       LINEAR_FIT_FUNCTION_3D_WITH_FIG_DIRECT(Pixel_units_fin,fig_switch);

    % Go beyond centroids that go beyond the edges
    % For simplicity make the terms greater than the ouside terms the outside term 

              r_3D(r_3D(:,1) > size_length(2),1) = size_length(2);
              r_3D(r_3D(:,2) > size_length(1),2) = size_length(1);

              xyq = r_3D;

           if fig_switch == 1
             figure, scatter3(Pixel_units_fin(:,1),Pixel_units_fin(:,2),(1:size(Pixel_units_fin,1))')
             hold on 
             scatter3(xyq(:,1),xyq(:,2),(1:size(xyq,1))','r'), grid on
           end

    else 
     % Here we need to make sure that all the sizes are the same        
            if size(Centroid_MID_cell_lump,2)==2
            xyq              =       [Centroid_MID_cell_lump, (1:size(Centroid_MID_cell_lump,1))'];      
            grad             =       min([180 - acosd(dot(normr_custom(xyq(1,:) - xyq(end,:)),[0 0 1]))  ;  acosd(dot(normr_custom(xyq(1,:) - xyq(end,:)),[0 0 1]))]); 

            else
            xyq  = Centroid_MID_cell_lump;
            grad = 0;

            end        

    end

    varargout{1} = xyq;
    varargout{2} = grad;


end
%% LINEAR_FIT_FUNCTION_3D_WITH_FIG_DIRECT
%   - Performs 3D linear fit on a set of points.
%   - Uses Cook's distance to exclude outliers.
%   - Returns fitted 3D coordinates and gradient along the principal direction.

function [varargout] = LINEAR_FIT_FUNCTION_3D_WITH_FIG_DIRECT(varargin)

Pixel_units_fin    =   varargin{1};
fig_switch         =   varargin{2};

% This function evaluates the turtuosity of the strands of the fibers and perfroms the split in terms of the gradient 
% of the sections of a fiber strand using a sectioning procedure analogous moving centroid approach 
% This function performs the following

%(a) extracts pixel locations which are outliers using the cooks distance 
%(b) extracts linearly independent points since svd works poorly on step slopes 
%(c) Performs svd on the points (under svd if steep slope occurs then maintain the points
    
% ====>>>>  This part of extracts only the pixels that pass the cook distance test
Centroid_MID_cell_lump              =     Pixel_units_fin;

if size(Centroid_MID_cell_lump,2)   == 2
 Centroid_MID_cell_lump_3D           =     [Centroid_MID_cell_lump, (1:size(Centroid_MID_cell_lump,1))'];
else
  Centroid_MID_cell_lump_3D           =    Centroid_MID_cell_lump;   
end

[centroid,potential_outlier]        =     regstats_func(Centroid_MID_cell_lump);

if ~isempty(centroid)
%  ====>>>>     Applying the number of linearly independent rows  to the variable centroid 
        tol = 1e-10;
        [~,rows_idx]                        =    LIROWS(centroid',tol);
        centroid(:,end+1)                   =    potential_outlier;
             
        if numel(rows_idx) > 1
% ====>>>>        Using the singular value decomposition approach 
         r0           =  mean(centroid);
         xyz_cent     =  bsxfun(@minus,centroid,r0);
         [~,~,V]      =  svd(xyz_cent,0);
         direction    =  V(:,1);
         direction    =  normc_custom(direction);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Step A : Straight forward approach
% ====>>>>       Correct for the straight lines 
                if any(direction == Inf | direction == -Inf | isnan(direction) | direction == -NaN )
                         r_3D                      =   [repmat(centroid(1,1:2), [size(Centroid_MID_cell_lump,1),1])  (1:size(Centroid_MID_cell_lump,1))'] ;
                         centroid_3D               =   Centroid_MID_cell_lump_3D;
                else                   
                         midpoints                 =   mean(Centroid_MID_cell_lump_3D);
                         dirvec                    =   direction';
                         
                         slice_num                 =   1 : size(Centroid_MID_cell_lump_3D,1);
                         dl                        =   (slice_num - midpoints(end))/dirvec(end);
                         
                         dl                        =    num2cell(dl,1);
                         points_slice              =    cellfun(@(x)([midpoints(1)  dirvec(1); midpoints(2)  dirvec(2)]*[1;x])',dl,'uni',0);
                         points_slice              =    vertcat(points_slice{:});
                         points_slice(:,end+1)     =    (1:size(points_slice,1))';
                         r_3D                      =    points_slice;                
                end 
% % Step B : will compute the perpendicular points on the line           
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ====>>>>      For pixels with little number independent rows  
        else 
                        centroid_3D     =   Centroid_MID_cell_lump_3D;
                        r_3D            =   Centroid_MID_cell_lump_3D;
        end

%   ====>>>>   This is due to the fact that the all the numbers are the same for the cook_distance thus the centroid is given as empty set   
else   
centroid_3D              =    Centroid_MID_cell_lump_3D;
r_3D                     =    Centroid_MID_cell_lump_3D;

end

% Cross check with the projection 
       if fig_switch == 1
          figure, scatter3(centroid_3D(:,1),centroid_3D(:,2),centroid_3D(:,3),'b','LineWidth',3)
          hold on 
          scatter3(r_3D(:,1),r_3D(:,2),r_3D(:,3),'r','filled'), grid on
          hold off
       end
      
       direct_vect      =       normr_custom(r_3D(1,:) - r_3D(end,:));
       grad             =       min([180 - acosd(dot(direct_vect,[0 0 1]))  ;  acosd(dot(direct_vect,[0 0 1]))]); 

       varargout{1} = r_3D;
       varargout{2} = grad;
end
%% STATS_BREAK_FUNCTION_22_update
%   - Processes metadata for each fiber segment.
%   - Extracts in-plane orientation and major/minor axes.
%   - Ensures consistency of fiber indexing after segmentation.

function [varargout] = STATS_BREAK_FUNCTION_22_update(varargin)
% This function breaks up the metadata of the STATS structure and extracts the orientation metadata from the variable 
% Also ensures that the index of the fibers are duly updated for the spatial instance assembler automator 

fin_lump_VAL_CELL_TEMP              =   varargin{1};
STATS_TEMP                          =   varargin{2};
GLOBAL_CENT_TOTAL_TEMP_refit        =   varargin{3};
GLOBAL_CENT_TOTAL_TEMP              =   varargin{4};
grad_out                            =   varargin{5};

fin_lump_VAL_CELL_temp_main     =    fin_lump_VAL_CELL_TEMP;
fields                          =    {'Centroid','Out_of_plane_Orientation','Final_Area','Pixel_List','Pixel_IDX_List'};

GLOBAL_CENT_TOTAL_TEMP_refit_UP   =   GLOBAL_CENT_TOTAL_TEMP_refit;
GLOBAL_CENT_TOTAL_TEMP_UP         =   GLOBAL_CENT_TOTAL_TEMP;
grad_out_UP                       =   grad_out;

% disp('initialize')
for ii = 1  :  length(fin_lump_VAL_CELL_temp_main)   
       if iscell(fin_lump_VAL_CELL_temp_main{ii})              
            if any(cell2mat(cellfun(@(x)iscell(x),fin_lump_VAL_CELL_temp_main{ii},'uni',0))==1)        
               fin_lump_VAL_CELL_temp_main{ii}      =   horzcat(fin_lump_VAL_CELL_temp_main{ii}{:});  
               GLOBAL_CENT_TOTAL_TEMP_refit_UP{ii}  =   horzcat(GLOBAL_CENT_TOTAL_TEMP_refit_UP{ii}{:});
               GLOBAL_CENT_TOTAL_TEMP_UP{ii}        =   horzcat(GLOBAL_CENT_TOTAL_TEMP_UP{ii}{:});
               grad_out_UP{ii}                      =   horzcat(grad_out_UP{ii}{:});
            end             
       end    
end
fin_lump_VAL_CELL_TEMP_O                         =   fin_lump_VAL_CELL_temp_main;
% disp('end initialize')

% Here assumptions to first and last are such that for fiber 1 and 2 we expect the cllusion to occur from the somewhere further from begining to end 
% Similarly assumptions to last and last but one are such that we expect the cllusion to occur from the somewhere from begining to ckoser to mid way


% disp('begin fin correction')

        In_plane_angle_CELL_FINAL                =         cell(1,length(fin_lump_VAL_CELL_temp_main)); 
        Maj_Min_CELL_FINAL                       =         cell(1,length(fin_lump_VAL_CELL_temp_main));

for ii = 1:length(fin_lump_VAL_CELL_temp_main)
    
    if iscell(fin_lump_VAL_CELL_temp_main{ii})
           
        for jj = 1:length(fin_lump_VAL_CELL_temp_main{ii})   
            if jj == 1
        f_idx                                            =    find(ismember(fin_lump_VAL_CELL_temp_main{ii}{1},fin_lump_VAL_CELL_temp_main{ii}{2}),1);       
        
        % though unlikely but to be safe  
        if isempty(f_idx)
        fin_lump_VAL_CELL_TEMP_O{ii}{1}                  =    fin_lump_VAL_CELL_temp_main{ii}{1};
        else 
        fin_lump_VAL_CELL_TEMP_O{ii}{1}(f_idx:end)       =    fin_lump_VAL_CELL_temp_main{ii}{1}(f_idx-1);   
        end 
               
        
        temp_STATS                                             =     STATS_TEMP{ii}(fin_lump_VAL_CELL_TEMP_O{ii}{1});
        temp_STATS                                             =     rmfield(temp_STATS,fields);        
        In_plane_angle                                         =     [temp_STATS.In_plane_Orientation];
        Maj_Min                                                =     [vertcat(temp_STATS.MajorAxisLength)  vertcat(temp_STATS.MinorAxisLength)];
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
            elseif jj == length(fin_lump_VAL_CELL_temp_main{ii})
                
        l_idx                                            =    find(ismember(fin_lump_VAL_CELL_temp_main{ii}{jj},fin_lump_VAL_CELL_temp_main{ii}{jj-1}),1,'last');
              
        if isempty(l_idx)
        fin_lump_VAL_CELL_TEMP_O{ii}{jj}                 =    fin_lump_VAL_CELL_temp_main{ii}{jj};
        else 
        fin_lump_VAL_CELL_TEMP_O{ii}{jj}(1:l_idx)        =    fin_lump_VAL_CELL_temp_main{ii}{jj}(l_idx+1);   
        end   
               
        temp_STATS                                             =    STATS_TEMP{ii}(fin_lump_VAL_CELL_TEMP_O{ii}{jj});  
        temp_STATS                                             =    rmfield(temp_STATS,fields);       
        In_plane_angle                                         =    [temp_STATS.In_plane_Orientation];
        Maj_Min                                                =    [vertcat(temp_STATS.MajorAxisLength)  vertcat(temp_STATS.MinorAxisLength)];
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
        
            elseif jj > 1 && jj < length(fin_lump_VAL_CELL_temp_main{ii})
                
        l_idx                                            =    find(ismember(fin_lump_VAL_CELL_temp_main{ii}{jj},fin_lump_VAL_CELL_temp_main{ii}{jj-1}),1,'last');
        f_idx                                            =    find(ismember(fin_lump_VAL_CELL_temp_main{ii}{jj},fin_lump_VAL_CELL_temp_main{ii}{jj+1}),1);
                    
       if isempty(l_idx) && any(f_idx ~= 1) 
            fin_lump_VAL_CELL_TEMP_O{ii}{jj}(f_idx:end)          =    fin_lump_VAL_CELL_temp_main{ii}{jj}(f_idx-1); 
       
        elseif isempty(f_idx) && any(l_idx ~= length(fin_lump_VAL_CELL_temp_main{ii}{jj}))            
            fin_lump_VAL_CELL_TEMP_O{ii}{jj}(1:l_idx)            =    fin_lump_VAL_CELL_temp_main{ii}{jj}(l_idx+1);
        
        elseif ~isempty(f_idx) && any(f_idx~=1)   && ~isempty(l_idx) && any(l_idx ~= length(fin_lump_VAL_CELL_temp_main{ii}{jj}))          
            fin_lump_VAL_CELL_TEMP_O{ii}{jj}(f_idx:end)          =    fin_lump_VAL_CELL_temp_main{ii}{jj}(f_idx-1);  
            fin_lump_VAL_CELL_TEMP_O{ii}{jj}(1:l_idx)            =    fin_lump_VAL_CELL_temp_main{ii}{jj}(l_idx+1);
       
        elseif  isempty(f_idx) && isempty(l_idx)       
            fin_lump_VAL_CELL_TEMP_O{ii}{jj}                     =    fin_lump_VAL_CELL_temp_main{ii}{jj};             
        
        elseif  any(f_idx ==1) || any(l_idx == 1)      
            fin_lump_VAL_CELL_TEMP_O{ii}{jj}                     =    fin_lump_VAL_CELL_temp_main{ii}{jj};  
            
        end
              
        temp_STATS                                             =    STATS_TEMP{ii}(fin_lump_VAL_CELL_TEMP_O{ii}{jj});    
        temp_STATS                                             =    rmfield(temp_STATS,fields);        
        In_plane_angle                                         =    [temp_STATS.In_plane_Orientation];
        Maj_Min                                                =    [vertcat(temp_STATS.MajorAxisLength)  vertcat(temp_STATS.MinorAxisLength)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                  
            end 
            

        In_plane_angle_CELL_FINAL{ii}{jj}              =         In_plane_angle; 
        Maj_Min_CELL_FINAL{ii}{jj}                     =         Maj_Min;
        end    
    else
        
        In_plane_angle_CELL_FINAL{ii}                  =         [STATS_TEMP{ii}.In_plane_Orientation];
        Maj_Min_CELL_FINAL{ii}                         =         [vertcat(STATS_TEMP{ii}.MajorAxisLength)  vertcat(STATS_TEMP{ii}.MinorAxisLength)];
    end
    
end

varargout{                  1}                  =    In_plane_angle_CELL_FINAL;
varargout{              end+1}                  =    fin_lump_VAL_CELL_temp_main;
varargout{              end+1}                  =    GLOBAL_CENT_TOTAL_TEMP_refit_UP;
varargout{              end+1}                  =    GLOBAL_CENT_TOTAL_TEMP_UP;
varargout{              end+1}                  =    grad_out_UP;
varargout{              end+1}                  =    Maj_Min_CELL_FINAL;

end
%% PIXEL_INTERPOLATE_FORMULATOR_SMART_ELLIPSE
%   - Populates pixel coordinates in ellipses around centroids.
%   - Adjusts ellipse sizes based on local fiber dimensions.
%   - Outputs pixel indices for each slice and checks for defective fibers.

function [varargout] = PIXEL_INTERPOLATE_FORMULATOR_SMART_ELLIPSE(varargin)

in_plane_angle            =    varargin{1};
centroid_track_temp       =    varargin{2};
size_length               =    varargin{3};
Maj_Min_CELL_value        =    varargin{4};

% In this step, we find the largest area within one standard of the mean we will use that to populate the ellipses on top and below
% Find the standard deviation between area values of reconstructed fiber 
% here make the adjusted diameter is minRadius which is the mean of values less than standard deviation away from the mean
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%

if numel(unique(Maj_Min_CELL_value(:,2))) == 1
   minRadius          =   mean(Maj_Min_CELL_value(:,2));
else 
   minRadius          =   mean(Maj_Min_CELL_value(Maj_Min_CELL_value(:,2) < mean(Maj_Min_CELL_value(:,2)) + std(Maj_Min_CELL_value(:,2)),2));    
end

if numel(unique(Maj_Min_CELL_value(:,1))) == 1
   majRadius          =   mean(Maj_Min_CELL_value(:,1));
else 
   majRadius          =   mean(Maj_Min_CELL_value(Maj_Min_CELL_value(:,1) < mean(Maj_Min_CELL_value(:,1)) + std(Maj_Min_CELL_value(:,1)),1));    
end

%%
% Note the minRadius and majRadius are slightly shrinked null any inducement of neighboring reigions during fiber stitching.
shrink_factor           =    0.90;
minRadius_temp          =    0.5*shrink_factor*minRadius;
majRadius_temp          =    0.5*shrink_factor*majRadius; 
Orientation_2D          =    deg2rad(in_plane_angle);

if size(Orientation_2D,1)==1
   Orientation_2D  = -Orientation_2D';
end

xCenter            =    num2cell(centroid_track_temp(:,1));
yCenter            =    num2cell(centroid_track_temp(:,2));
phie               =    num2cell(Orientation_2D);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Constrain the values between (but including) zero and the max size of the matrix by one 
% c = cellfun(@(x)x(:,x(1,:)>0 & x(1,:)<size_length_2D(2)+1),c,'uni',0);   % c has numbers iin ncreasing order along each cols 
% r = cellfun(@(y)y(y(:,1)>0 & y(:,1)<size_length_2D(1)+1,:),r,'uni',0);   % r has numbers in increasing order along each rows

[c,r]= cellfun(@(x,y)meshgrid(round(x-majRadius_temp):round(x+majRadius_temp),round(y-majRadius_temp):round(y+majRadius_temp)),xCenter,yCenter,'uni',0);

% Pseudo_Linear_Index 
Pixel_Idx = cellfun(@(x,y,z,c,r) find( ((c - x)*cos(z) + (r - y)*sin(z)).^2 /(majRadius_temp)^2  ...
                   +(-(c - x)*sin(z) + (r - y)*cos(z)).^2 /(minRadius_temp)^2  <= 1) , ...
                       xCenter, yCenter, phie,c,r,'uni',0);
                  
% Actual Pixel_List
Pixel_List = cellfun(@(x,r,c) [r(x) c(x)], Pixel_Idx, r, c,'uni',0); % Here its in rows and columns 
Pixel_List = cellfun(@(x)fliplr(x),Pixel_List,'uni',0);

% % Replace the empty cells with the original Pixel values 
empty_idx = find(cell2mat(cellfun(@(x)isempty(x),Pixel_List,'uni',0)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if find(ismember(1:length(Pixel_Idx),empty_idx),1,'first') == 1 
%    [Pixel_List{1:find(~ismember(1:length(Pixel_Idx),empty_idx),1,'first')}] = deal(Pixel_List{find(~ismember(1:length(Pixel_Idx),empty_idx),1,'first')});
%    
%    
% elseif  find(ismember(1:length(Pixel_Idx),empty_idx),1,'last') == length(Pixel_Idx)
%    [Pixel_List{find(~ismember(1:length(Pixel_Idx),empty_idx),1,'last'):end}] = deal(Pixel_List{find(~ismember(1:length(Pixel_Idx),empty_idx),1,'last')});
% 
% else 
%    [Pixel_List{1:find(~ismember(1:length(Pixel_Idx),empty_idx),1,'first')}]  = deal(Pixel_List{find(~ismember(1:length(Pixel_Idx),empty_idx),1,'first')});
%    [Pixel_List{find(~ismember(1:length(Pixel_Idx),empty_idx),1,'last'):end}] = deal(Pixel_List{find(~ismember(1:length(Pixel_Idx),empty_idx),1,'last')});
% end

% Eliminate defective fibers 
if isempty(empty_idx)
   defective = 0;
    % Contrain the Pixel list components 
    Pixel_List = cellfun(@(x)x((ismember(x(:,1),1:size_length(2)) & ismember(x(:,2),1:size_length(1))),:),Pixel_List,'uni',0);

    % Actual Pixel_IDX_List
    Pixel_Idx = cellfun(@(x)sub2ind([size_length(1)   size_length(2)],x(:,2),x(:,1)),Pixel_List,'uni',0);


    for ii = 1:length(Pixel_Idx)
    [STATS_NEW(ii).FINAL_VALUES]         =   Pixel_List{ii};
    [Pixel_Idx_FINAL_VALUES{ii}]         =   Pixel_Idx{ii};
    end
    
    minRadius_idxval                  =     [minRadius length(STATS_NEW)];

else
    defective = 1;
    Pixel_Idx_FINAL_VALUES    =   [];
    STATS_NEW                 =   [];
    minRadius_idxval          =   [];
       
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


varargout{               1}       =     Pixel_Idx_FINAL_VALUES;
varargout{           end+1}       =     STATS_NEW;
varargout{           end+1}       =     minRadius_idxval;
varargout{           end+1}       =     defective;
end


%% SPATIAL_INSTANCE_ASSEMBLER_AUTOMATOR_MULTI
%   - Converts 2D pixel coordinates into 3D linear indices in the volume.
%   - Stores the corresponding height (slice) for each pixel.
%   - Handles both single and multi-fiber cases.

function [varargout] = SPATIAL_INSTANCE_ASSEMBLER_AUTOMATOR_MULTI(varargin)

STATS_NEW_MIDDLE_LINE            =   varargin{1};
fin_lump_VAL_CELL_TEMP_O         =   varargin{2};
iter_slice                       =   varargin{3};
size_length                      =   varargin{4};

if iscell(STATS_NEW_MIDDLE_LINE)
     check_mat_center       =   cell(1,length(STATS_NEW_MIDDLE_LINE));
     height_variable        =   cell(1,length(STATS_NEW_MIDDLE_LINE));
     
     for ii = 1:length(STATS_NEW_MIDDLE_LINE)        % Disregarding the defective cell;
         STATS_NEW_MIDDLE_LINE_temp    =    STATS_NEW_MIDDLE_LINE{ii};
         
         if ~isempty(STATS_NEW_MIDDLE_LINE_temp)
         hh                            =    iter_slice + fin_lump_VAL_CELL_TEMP_O{ii}(1)-1;

             Height_Info_slices  = cell(1,length(STATS_NEW_MIDDLE_LINE_temp));
             for kk = 1:length(STATS_NEW_MIDDLE_LINE_temp)
             [STATS_NEW_MIDDLE_LINE_temp(kk).Pixel_IDX_List]   =  deal (sub2ind ([size_length(1) size_length(2) size_length(3)],STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES(:,end),...
                                                                                  STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES(:,1),(hh-1+kk)*ones(size(STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES,1),1)));
             
             
             Height_Info_slices{kk} = (hh-1+kk)*ones(size(STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES,1),1);                                                                 
                                                                              
             end

         check_mat_center{ii}                           =     vertcat(STATS_NEW_MIDDLE_LINE_temp.Pixel_IDX_List);
         height_variable{ii}                            =     vertcat(Height_Info_slices{:});
         end

     end

else 
        STATS_NEW_MIDDLE_LINE_temp   =   STATS_NEW_MIDDLE_LINE;
        
        if ~isempty(STATS_NEW_MIDDLE_LINE_temp)
        hh                           =   iter_slice;

        Height_Info_slices  = cell(1,length(STATS_NEW_MIDDLE_LINE_temp));
        for kk = 1:size(STATS_NEW_MIDDLE_LINE_temp,2)
        [STATS_NEW_MIDDLE_LINE_temp(kk).Pixel_IDX_List]   =  deal (sub2ind ([size_length(1) size_length(2) size_length(3)],STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES(:,end),...
                                                                            STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES(:,1),(hh-1+kk)*ones(size(STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES,1),1)));
                                                                        
         Height_Info_slices{kk} = (hh-1+kk)*ones(size(STATS_NEW_MIDDLE_LINE_temp(kk).FINAL_VALUES,1),1);                                                                
                                                                        
        end

        check_mat_center                           =   vertcat(STATS_NEW_MIDDLE_LINE_temp.Pixel_IDX_List);
        height_variable                            =   vertcat(Height_Info_slices{:});
        else 
        check_mat_center                           =   [];
        height_variable                            =   [];               
        end
end

varargout{      1} = check_mat_center;
varargout{  end+1} = height_variable;

end 

%% regstats_func
%   - Identifies outlier points using regression statistics and Cook's distance.
%   - Returns the filtered set of 3D coordinates for linear fitting.

function [varargout]  = regstats_func(varargin)

Centroid_MID_cell_lump = varargin{1};
xy  =  fliplr(Centroid_MID_cell_lump);

% figure, scatter(xy(:,end),xy(:,1),'*')
% stats = regstats(Y,X,'linear');

stats = regstats(xy(:,end),xy(:,1),'linear');
% if Cook's Distance > n/4 is a typical treshold that is used to suggest the presence of an outlier

% potential_outlier = find(mean(stats.cookd)-0.4*std(stats.cookd) < stats.cookd  & stats.cookd < mean(stats.cookd) + 0.4*std(stats.cookd));
potential_outlier = find(mean(stats.cookd)-std(stats.cookd) < stats.cookd  & stats.cookd < mean(stats.cookd) + std(stats.cookd));

fin_out = fliplr(xy(potential_outlier,:));
varargout{1} = fin_out;
varargout{2} = potential_outlier;
end

%% LIROWS
%   - Extracts linearly independent rows from a matrix.
%   - Used to ensure stable SVD computations during linear fitting.

function [Xsub,idx]=LIROWS(X,tol)
%Extract a linearly independent set of columns of a given matrix X
%    [Xsub,idx]=licols(X)
%in:
%  X: The given input matrix
%  tol: A rank estimation tolerance. Default=1e-10
%out:
% Xsub: The extracted columns of X
% idx:  The indices (into X) of the extracted columns
     if ~nnz(X) %X has no non-zeros and hence no independent columns
         Xsub=[]; idx=[];
         return
     end
     if nargin<2, tol=1e-10; end
       [Q, R, E] = qr(X,0); 
       if ~isvector(R)
        diagr = abs(diag(R));
       else
        diagr = R(1);   
       end
       %Rank estimation
       r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
       idx=sort(E(1:r));
       Xsub=X(:,idx);
       Xsub = Xsub';
end

%% RETRACKER_FUNCTION
%   - Maps fiber slices to their correct positions in the global volume.
%   - Returns starting index and slice range for assembly.

function [ii_idx , slice_region] = RETRACKER_FUNCTION(numel_ellipses,length_check)

     for ii = 1:length(numel_ellipses)
         slice_region = find(numel_ellipses{1,ii}==length_check);
         if ~isempty(slice_region)
         break
         end
     end
     
ii_idx = numel_ellipses{2,ii};

end
