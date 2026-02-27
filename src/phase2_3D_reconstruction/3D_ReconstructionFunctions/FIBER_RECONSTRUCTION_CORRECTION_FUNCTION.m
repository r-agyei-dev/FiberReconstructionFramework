%% This latest version accounts for defective fibers removed from the dataset

% This function performs: 
% (a) Separation of fibers using a dedicated algorithm
% (b) Radial population of fibers based on the centroid pixel locations
% (c) Verification to ensure all segmented areas are removed before fiber stitching
% (d) Creation of an HDF file for fibers prior to stitching
% (e) 3D stitching of the fiber volume
% (f) Interpolation to refine the reconstructed 3D fiber volume

% Functions utilized include: 
%  FIBER_SEPERATING_ALGORITHM_CODE_DENSE_FMT_UPDATED
%  CENTRAL_PIXEL_GLOBAL_ESTIMATOR_UPDATED_DIRECT_DENSE_v2_TT_EX
%  SPLITTER_CORRECTOR_FUNC_A
%  Splitter_Lumper_func
%  THREE_D_LABELLER_HDF_CREATOR_TT
%  THREE_D_STITCHER_HYBRID_v2_CONDENSED_FORMAT_TT
%  PIXEL_INTERPOLATE_FORMULATOR_SMART_ELLIPSE_FURTHERv2_TT_EX

function [varargout] = FIBER_RECONSTRUCTION_CORRECTION_FUNCTION(varargin)

% Input arguments:
% properties                 = Cell array containing metadata of ellipses
% iter_limit                 = Maximum number of slices to stack ellipses
% size_length                = Size of the 3D reconstruction volume
% load_variable              = Indicates the axis of the serial sections for the 3D volume
% angular_thresh             = Threshold for out-of-plane fiber angle
% save_path                  = Directory to save 3D reconstructed volumes
% tomo_dataset_idx           = Index of the tomographic dataset
% save_path_extra_file       = (Optional) Extra save directory for reconstructed volumes

properties                      =   varargin{1};
iter_limit                      =   varargin{2};
size_length                     =   varargin{3};
load_variable                   =   varargin{4};
angular_thresh                  =   varargin{5};
save_path                       =   varargin{6};
tomo_dataset_idx                =   varargin{7};
% save_path_extra_file          =   varargin{8};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1: Label each fiber and associate with its slice
disp('==> STEP 1: Label fibers and associate with slice index')
% fib_idx_n_slice_loc_cell stores: 
%   Row 1: consecutive fiber indices
%   Row 2: corresponding slice numbers

empty_cell_idx                  = find(cell2mat(cellfun(@(x)~isempty(x),properties,'uni',0)));    
fib_idx_n_slice_loc_cell        = cell(2,iter_limit); 

for aa = 1:length(empty_cell_idx)
    if aa == 1
        fib_idx_n_slice_loc_cell{1,aa} = 1:length(properties{empty_cell_idx(aa)});   
    else
        fib_idx_n_slice_loc_cell{1,aa} = 1 + fib_idx_n_slice_loc_cell{1,(aa-1)}(end) : ...
                                        fib_idx_n_slice_loc_cell{1,(aa-1)}(end)+length(properties{empty_cell_idx(aa)});
    end        
    fib_idx_n_slice_loc_cell{2,aa} = empty_cell_idx(aa);        
end  

clearvars SLICE_REGIONS_STATS_STRUCTURE Centroid_cell track_MAIN_CELL sieve_parameter_unique_CELL    
clearvars HDF_ARRAY_MAT temp_HDF_MD_ARRAY 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2: Remove tiny segments below threshold
disp('==> STEP 2: Remove tiny segments')
thresh_for_noise         = 2;
main_cell                = [properties{:}]; % Flatten all properties
length_cell              = cell2mat(cellfun(@(x)length(x),main_cell,'uni',0)); % Fiber lengths
index_to_check           = find(length_cell > thresh_for_noise); % Keep fibers above threshold
index_to_check(2,:)      = length_cell(index_to_check); % Store lengths of surviving fibers
idx_after_first_thresh   = index_to_check(1,:);
main_cell_update         = main_cell(idx_after_first_thresh);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 3: Execute fiber separating algorithm
disp('==> STEP 3: Execute fiber separation')
% Outputs:
% Distinct_Fiber_IDX_CELL         = indices of consecutive ellipses per fiber
% GLOBAL_CENT_Distinct_Fiber      = centroids of consecutive ellipses
% ROGUE_EXTR, ROGUE_BEGIN_SLICE   = info about rogue regions

[Distinct_Fiber_IDX_CELL,GLOBAL_CENT_Distinct_Fiber,~,~,~,~] = ...
    FIBER_SEPERATING_ALGORITHM_FUNCTION(idx_after_first_thresh,fib_idx_n_slice_loc_cell,properties,size_length);      

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 4: Radially populate fibers at centroid pixel locations
disp('==> STEP 4: Populate fibers based on centroid pixels')
% Outputs:
% Linear_Index_center_GLOBAL   = logical indices of reconstructed fibers
% ORIGI_Centroid_MID_cel       = original centroid coordinates
% Pixel_Idx_cell_GLOBAL        = per-slice pixel indices
% Sub_Idx_Split_Volume         = split fiber indices
% height_variable              = slice indices

temp_thresh = 1:length(idx_after_first_thresh);
Distinct_Fiber_IDX_CELL_temp = Distinct_Fiber_IDX_CELL(temp_thresh);
GLOBAL_CENT_Distinct_Fiber_temp = GLOBAL_CENT_Distinct_Fiber(temp_thresh);
idx_after_first_slice_temp = idx_after_first_thresh(temp_thresh);

tic
[~,Linear_Index_center_GLOBAL,~,ORIGI_Centroid_MID_cel,...
 Pixel_Idx_cell_GLOBAL,~,Sub_Idx_Split_Volume,~,height_variable] = ...
    CENTRAL_PIXEL_POPULATOR_FUNCTION(Distinct_Fiber_IDX_CELL_temp,GLOBAL_CENT_Distinct_Fiber_temp,...
                                     main_cell,fib_idx_n_slice_loc_cell,idx_after_first_slice_temp,...
                                     size_length,temp_thresh);                                                                                                                                                                                                                                               
toc

disp('Remove defective fibers')
% Remove entries corresponding to empty fibers
defective_idx = cell2mat(cellfun(@(x)isempty(x),Linear_Index_center_GLOBAL,'uni',0));
Linear_Index_center_GLOBAL(defective_idx == 1) = []; 
ORIGI_Centroid_MID_cel(defective_idx == 1) = []; 
Pixel_Idx_cell_GLOBAL(defective_idx == 1) = []; 
Sub_Idx_Split_Volume(defective_idx == 1) = []; 
height_variable(defective_idx == 1) = []; 
main_cell_update(defective_idx == 1) = []; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 5: Update fiber properties based on splits
disp('Update fiber properties after splits')
In_plane_extract          = cell(1,size(Sub_Idx_Split_Volume,2));
major_axis_length_extract = cell(1,size(Sub_Idx_Split_Volume,2));
minor_axis_length_extract = cell(1,size(Sub_Idx_Split_Volume,2));

for ii = 1:size(Sub_Idx_Split_Volume,2)
    if ~iscell(Sub_Idx_Split_Volume{ii})
        In_plane_extract{ii} = [main_cell_update{ii}.In_plane_Orientation]';
        major_axis_length_extract{ii} = [main_cell_update{ii}.MajorAxisLength]';
        minor_axis_length_extract{ii} = [main_cell_update{ii}.MinorAxisLength]';
    else
        for tt = 1:size(Sub_Idx_Split_Volume{ii},2)
            In_plane_extract{ii}{tt} = [main_cell_update{ii}(Sub_Idx_Split_Volume{ii}{tt}).In_plane_Orientation]';  
            major_axis_length_extract{ii}{tt} = [main_cell_update{ii}(Sub_Idx_Split_Volume{ii}{tt}).MajorAxisLength]';
            minor_axis_length_extract{ii}{tt} = [main_cell_update{ii}(Sub_Idx_Split_Volume{ii}{tt}).MinorAxisLength]';
        end
    end
end
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 5A: Ensure all segmented areas removed before stitching
disp('===>> STEP 5: Remove all segmented regions before stitching')
[Linear_Index_center_GLOBAL,minor_axis_length_extract,...
 ORIGI_Centroid_MID_cel,Pixel_Idx_cell_GLOBAL,numel_count] = ...
    SPLITTER_CORRECTOR_FUNC_A(Linear_Index_center_GLOBAL,minor_axis_length_extract,ORIGI_Centroid_MID_cel,...
                               height_variable,Pixel_Idx_cell_GLOBAL);                                                                

% Remove empty cells from all variables
empty_cells = cell2mat(cellfun(@(x)isempty(x),Linear_Index_center_GLOBAL,'uni',0));
Linear_Index_center_GLOBAL(empty_cells) = [];
numel_count(empty_cells) = [];
ORIGI_Centroid_MID_cel(empty_cells) = [];
minor_axis_length_extract(empty_cells) = [];
Pixel_Idx_cell_GLOBAL(empty_cells) = [];

% Lump fibers and centroids for HDF creation
[Linear_Index_center_GLOBAL_LUMP] = Splitter_Lumper_func(Linear_Index_center_GLOBAL);                                                           
[ORIGI_Centroid_MID_cel_LUMP] = Splitter_Lumper_func(ORIGI_Centroid_MID_cel);    
[minor_axis_length_extract_LUMP] = Splitter_Lumper_func(minor_axis_length_extract);  
Centroid_MID_cell_LUMP_ORIG = ORIGI_Centroid_MID_cel_LUMP;

Centroid_MID_cell_LUMP_ORIG(cell2mat(cellfun(@(x)isempty(x),Centroid_MID_cell_LUMP_ORIG,'uni',0))) = [];
minor_axis_length_extract_LUMP(cell2mat(cellfun(@(x)isempty(x),minor_axis_length_extract_LUMP,'uni',0))) = [];
Linear_Index_center_GLOBAL_LUMP(cell2mat(cellfun(@(x)isempty(x),Linear_Index_center_GLOBAL_LUMP,'uni',0))) = [];
[Pixel_Idx_cell_GLOBAL_LUMP] = Pixel_Idx_cell_GLOBAL_LUMPER(Pixel_Idx_cell_GLOBAL,numel_count);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 6: Create HDF file for fibers without stitching
disp('==> STEP 6: Create HDF file for fibers before stitching')
lump = 0;       % Idx_height not required; volume labeled as is
hdfsave = 1;
flip_query = 1;

[TEMP_MAT,~] = THREE_D_LABELLER_HDF_CREATOR_A(Linear_Index_center_GLOBAL_LUMP,size_length,lump,...
                                               hdfsave,Centroid_MID_cell_LUMP_ORIG,flip_query,...
                                               load_variable,save_path,tomo_dataset_idx); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 7: Filter fibers by orientation
disp('===>> STEP 7: Filter fibers by out-of-plane orientation')
% Compute in-plane angles and out-of-plane angles to identify good and bad fibers
Exclude_bad = 0;
in_plane_cell = cellfun(@(x) atand((x(1,2)-x(end,2))/(x(end,1)-x(1,1))),Centroid_MID_cell_LUMP_ORIG ,'uni',0); 
single_temp = cell2mat(cellfun(@(x)size(x,1),Centroid_MID_cell_LUMP_ORIG ,'uni',0)) == 1;
[in_plane_cell{single_temp}] = deal(0);

grad_out_LUMP = cell2mat(cellfun(@(x,p) 90 - acosd((x(end,1)-x(1,1))/(norm(x(end,:)-x(1,:))*cosd(p))),...
                                  Centroid_MID_cell_LUMP_ORIG ,in_plane_cell,'uni',0));
grad_out_LUMP = num2cell(grad_out_LUMP);
[grad_out_LUMP{cell2mat(cellfun(@(x)isnan(x),grad_out_LUMP,'uni',0))}] = deal(0);
grad_out_LUMP = cell2mat(grad_out_LUMP);

good_fibers = find(abs(grad_out_LUMP)<angular_thresh);
bad_fibers = find(~ismember(1:length(grad_out_LUMP),good_fibers));

Exclude_bad = 0;
if Exclude_bad == 1
    TEMP_MAT(~ismember(TEMP_MAT,good_fibers)) = 0;
else 
    positive_nonzeros = unique(nonzeros(TEMP_MAT));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 8: Stitch over-segmented fibers
disp('===>> STEP 8: Stitch over-segmented fibers')
% Determine which fibers to stitch based on filtering
if Exclude_bad == 1
    avail_terms = good_fibers;
else 
    avail_terms = positive_nonzeros;
end

ang_thresh = 20;       % maximum out-of-plane angle allowed
inplane_thresh = 5;    % maximum in-plane deviation allowed

[IDX_to_keep,fib_diam_stitched,Linear_Index_stitched,...
 FINAL_CENTROID_stitched,min_slice_stitched] = ...
    THREE_D_STITCHER_OF_SEG_FIBERS_FUNCTION(avail_terms,Linear_Index_center_GLOBAL_LUMP,...
                                            size_length,minor_axis_length_extract_LUMP,...
                                            Centroid_MID_cell_LUMP_ORIG);                                                                                                        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 9: Pixel population after fiber stitching
disp('===>> STEP 9: Populate pixels for stitched fibers')
multi_indices  = find(cell2mat(cellfun(@(x)any(numel(x)>1),IDX_to_keep,'uni',0)));
Pixel_Idx_cell_GLOBAL_LUMP_PART = Pixel_Idx_cell_GLOBAL_LUMP(cell2mat(IDX_to_keep(cell2mat(cellfun(@(x)numel(x)==1,IDX_to_keep,'uni',0)))));

tic
[Pixel_Idx_Per_Slice_MAIN,height_variable_MAIN,~,...
 Linear_Index_Final,OUT_PLANE_ANGLE,FINAL_CORRECTED_CENTROID,...
 fib_diam_stitched_MAIN_VARIABLE] = ...
    PIXEL_INTERPOLATE_FORMULATOR_SMART_ELLIPSE_FURTHER_B(fib_diam_stitched,FINAL_CENTROID_stitched,...
                                                          size_length,Linear_Index_stitched,min_slice_stitched,...
                                                          multi_indices,Pixel_Idx_cell_GLOBAL_LUMP_PART);                                                                                                                                                                                                                                                                                                                                                                                                                   
toc 

% Return outputs
varargout{1} = Linear_Index_Final;
varargout{end+1} = Pixel_Idx_Per_Slice_MAIN;
varargout{end+1} = FINAL_CORRECTED_CENTROID;
varargout{end+1} = IDX_to_keep;
varargout{end+1} = OUT_PLANE_ANGLE;
varargout{end+1} = height_variable_MAIN;
varargout{end+1} = fib_diam_stitched_MAIN_VARIABLE;