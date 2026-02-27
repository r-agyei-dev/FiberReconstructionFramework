% PIXEL_INTERPOLATE_FORMULATOR_SMART_ELLIPSE_FURTHER_B
% -------------------------------------------------------------------------
% Refines fiber reconstructions by fitting oriented ellipses along each
% centroid track and generating the corresponding voxel indices.
%
% Purpose:
%   • Uses stitched fiber centroids and diameters
%   • Estimates fiber orientation in 3D
%   • Expands each centroid into an ellipse cross-section
%   • Generates voxel indices representing the fiber volume
%   • Removes defective and overwritten fibers
%
% High-level workflow:
%   1) Select fibers that were stitched across slices
%   2) Compute in-plane and out-of-plane orientation angles
%   3) Convert diameters → ellipse major/minor radii
%   4) Rasterize ellipses at each centroid location
%   5) Assemble full 3D fiber indices
%   6) Merge with unstitched fibers
%   7) Remove overlaps (overwrite cleanup)
%
% Inputs (via varargin):
%   {1} fib_diam_stitched               - Estimated fiber diameters
%   {2} FINAL_CENTROID_stitched         - Cell array of centroid tracks
%   {3} size_length                     - Volume size [rows cols slices]
%   {4} Linear_Index_stitched           - Existing linear indices
%   {5} min_slice_stitched              - Starting slice per fiber
%   {6} multi_indices                   - Indices of stitched fibers
%   {7} Pixel_Idx_cell_GLOBAL_LUMP_PART - Unstitched pixel index cells
%
% Outputs (varargout):
%   1) Pixel indices per slice per fiber
%   2) Height (slice) coordinates
%   3) STATS_NEW_CELL (per-slice stats)
%   4) Final linear indices per fiber
%   5) Out-of-plane angles
%   6) Final centroid tracks
%   7) Fiber diameters

function [varargout] = PIXEL_INTERPOLATE_FORMULATOR_SMART_ELLIPSE_FURTHER_B(varargin)

fib_diam_stitched                       =      varargin{1};
FINAL_CENTROID_stitched                 =      varargin{2};
size_length                             =      varargin{3};
Linear_Index_stitched                   =      varargin{4};
min_slice_stitched                      =      varargin{5};
multi_indices                           =      varargin{6};
Pixel_Idx_cell_GLOBAL_LUMP_PART         =      varargin{7};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Select only the fibers that were stitched across slices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

FINAL_CENTROID_stitched_PART            =      FINAL_CENTROID_stitched(multi_indices);   
fib_diam_stitched_PART                  =      fib_diam_stitched(multi_indices);
Linear_Index_stitched_PART              =      Linear_Index_stitched(multi_indices);
min_slice_stitched_PART                 =      min_slice_stitched(multi_indices);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Compute orientation of each fiber from its endpoints
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Extract first/last centroid of each fiber
endpoints     =   cellfun(@(x)[x(1,:) ; x(end,:)],FINAL_CENTROID_stitched,'uni',0);

% Identify single-slice fibers
single_temp   =   cell2mat(cellfun(@(x)size(x,1),FINAL_CENTROID_stitched,'uni',0)) == 1;

% In-plane orientation angle (XY projection)
in_plane_cell = cellfun(@(x) ...
    atand((x(1,2)-x(end,2))/(x(end,1)-x(1,1))), ...
    FINAL_CENTROID_stitched,'uni',0);

% Force sensible defaults for degenerate cases
[in_plane_cell{single_temp}] = deal(0);
[in_plane_cell{cell2mat(cellfun(@(x) ...
    norm(x(1,1:2)-x(end,1:2)),FINAL_CENTROID_stitched,'uni',0)) == 0}] = deal(0);

% Out-of-plane tilt angle (3D inclination)
out_plane_angle = cellfun(@(x,p) ...
    90 - acosd((x(end,1)-x(1,1)) / ...
    (norm(x(end,:)-x(1,:))*cosd(p))), ...
    FINAL_CENTROID_stitched,in_plane_cell,'uni',0);

% Replace NaNs with zero tilt
[out_plane_angle{cell2mat(cellfun(@(x)isnan(x),out_plane_angle,'uni',0))}] = deal(0);
out_plane_angle = cell2mat(out_plane_angle);

Orientation_2D = cell2mat(in_plane_cell);
if size(Orientation_2D,1)==1
    Orientation_2D  = Orientation_2D';
end

out_plane_angle_PART =  out_plane_angle(multi_indices);   

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Convert diameter → ellipse radii (accounting for tilt)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

majRadius_vec = fib_diam_stitched_PART ./ cosd(out_plane_angle_PART);
minRadius_vec = fib_diam_stitched_PART;

Pixel_Idx_Per_Slice_MAIN = cell(1,length(FINAL_CENTROID_stitched_PART));
STATS_NEW_CELL           = cell(1,length(FINAL_CENTROID_stitched_PART)); 
defective                = zeros(1,length(FINAL_CENTROID_stitched_PART));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Rasterize ellipses along each centroid track
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for kk = 1:length(FINAL_CENTROID_stitched_PART)
    
    centroid_track_temp = FINAL_CENTROID_stitched_PART{kk};
    xCenter             = num2cell(centroid_track_temp(:,1));
    yCenter             = num2cell(centroid_track_temp(:,2));
    phie                = 180 - Orientation_2D(kk);
    minRadius           = 0.5 * minRadius_vec(kk);
    majRadius           = 0.5 * majRadius_vec(kk);                                                
                                                         
    % Build local grids around each centroid
    [c,r]= cellfun(@(x,y) ...
        meshgrid(round(x-majRadius):round(x+majRadius), ...
                 round(y-majRadius):round(y+majRadius)), ...
        xCenter,yCenter,'uni',0);

    % Ellipse membership test in rotated frame
    Pixel_Idx = cellfun(@(x,y,c,r) ...
        find(((c - x)*cos(phie) + (r - y)*sin(phie)).^2 /(majRadius)^2 + ...
             (-(c - x)*sin(phie) + (r - y)*cos(phie)).^2 /(minRadius)^2 <= 1), ...
        xCenter, yCenter,c,r,'uni',0);                                                      
                                                          
    % Convert to (col,row) pixel coordinates
    Pixel_List = cellfun(@(x,r,c) [r(x) c(x)], Pixel_Idx, r, c,'uni',0);
    Pixel_List = cellfun(@(x)fliplr(x),Pixel_List,'uni',0);                                                      
                                                          
    empty_idx = find(cell2mat(cellfun(@(x)isempty(x),Pixel_List,'uni',0)));                                                 

    % --- Reject defective fibers
    if isempty(empty_idx)

        % Clip to image bounds
        Pixel_List = cellfun(@(x) ...
            x((ismember(x(:,1),1:size_length(2)) & ...
               ismember(x(:,2),1:size_length(1))),:), ...
            Pixel_List,'uni',0);

        % Convert to linear indices
        Pixel_Idx = cellfun(@(x) ...
            sub2ind([size_length(1) size_length(2)],x(:,2),x(:,1)), ...
            Pixel_List,'uni',0);

        for ii = 1:length(Pixel_Idx)
            STATS_NEW(ii).FINAL_VALUES = Pixel_List{ii};
            Pixel_Idx_FINAL_VALUES{ii} = Pixel_Idx{ii};
        end                                                      

        Pixel_Idx_Per_Slice_MAIN{kk} = Pixel_Idx_FINAL_VALUES;
        STATS_NEW_CELL{kk}           = STATS_NEW;  

    else 
        defective(kk) = 1;
    end                                                        

    clearvars STATS_NEW Pixel_Idx_FINAL_VALUES
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Assemble full 3D fiber indices
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Linear_Index_Final = cell(1,length(Linear_Index_stitched_PART));
for kk = 1:length(Linear_Index_stitched_PART)   
    if defective(kk) == 0
        Linear_Index_Final{kk} = ...
            SPATIAL_INSTANCE_ASSEMBLER_AUTOMATOR_MULTI_FURTHER( ...
            STATS_NEW_CELL{kk}, ...
            min_slice_stitched_PART(kk), ...
            size_length);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Remove empty / defective fibers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Defective_fiber_Idx = cell2mat(cellfun(@(x)isempty(x),Linear_Index_Final,'uni',0));
Pixel_Idx_Per_Slice_MAIN(Defective_fiber_Idx)      = [];
STATS_NEW_CELL(Defective_fiber_Idx)                = [];
Linear_Index_Final(Defective_fiber_Idx)            = [];
out_plane_angle_PART(Defective_fiber_Idx)          = [];
FINAL_CENTROID_stitched_PART(Defective_fiber_Idx)  = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Merge back unstitched fibers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fib_diam_stitched(multi_indices)          = [];
FINAL_CENTROID_stitched(multi_indices)    = [];
Linear_Index_stitched(multi_indices)      = [];
out_plane_angle(multi_indices)            = [];

Linear_Index_Final_MAIN_VARIABLE           = [Linear_Index_stitched Linear_Index_Final];
out_plane_angle_MAIN_VARIABLE              = [out_plane_angle out_plane_angle_PART];
FINAL_CENTROID_stitched_PART_MAIN_VARIABLE = [FINAL_CENTROID_stitched FINAL_CENTROID_stitched_PART];
Pixel_Idx_Per_Slice_MAIN_VARIABLE          = [Pixel_Idx_cell_GLOBAL_LUMP_PART Pixel_Idx_Per_Slice_MAIN];
fib_diam_stitched_MAIN_VARIABLE            = [fib_diam_stitched fib_diam_stitched_PART];

disp('Creating the height variable for the respective fibers')
[~,~,height_variable_MAIN_VARIABLE] = ...
    cellfun(@(x)ind2sub(size_length,x), ...
    Linear_Index_Final_MAIN_VARIABLE,'uni',0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Remove overwritten fibers in the reconstructed volume
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

IMAGE_VOL_1 = zeros(size_length); 

disp('Begin Ridding all the overwritten fibers')
for ii = 1:length(Linear_Index_Final_MAIN_VARIABLE)
    IMAGE_VOL_1(Linear_Index_Final_MAIN_VARIABLE{ii}) = ii;
end
disp('End Ridding all the overwritten fibers')

disp('Remove all the overwritten fibers')
empty_indexes  = ~ismember( ...
    1:length(Linear_Index_Final_MAIN_VARIABLE), ...
    unique(nonzeros(IMAGE_VOL_1)));

Pixel_Idx_Per_Slice_MAIN_VARIABLE(empty_indexes)      = [];
height_variable_MAIN_VARIABLE(empty_indexes)          = [];
Linear_Index_Final_MAIN_VARIABLE(empty_indexes)       = [];
out_plane_angle_MAIN_VARIABLE(empty_indexes)          = [];
FINAL_CENTROID_stitched_PART_MAIN_VARIABLE(empty_indexes) = [];
fib_diam_stitched_MAIN_VARIABLE(empty_indexes)        = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% --- Outputs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout{1}       = Pixel_Idx_Per_Slice_MAIN_VARIABLE;
varargout{end+1}   = height_variable_MAIN_VARIABLE;
varargout{end+1}   = STATS_NEW_CELL;
varargout{end+1}   = Linear_Index_Final_MAIN_VARIABLE;
varargout{end+1}   = out_plane_angle_MAIN_VARIABLE;
varargout{end+1}   = FINAL_CENTROID_stitched_PART_MAIN_VARIABLE;
varargout{end+1}   = fib_diam_stitched_MAIN_VARIABLE;

end