%--------------------------------------------------------------------------
% FINAL_STATS_FORMULATOR_SMART_INTERSECT_CORRECTOR
%
% PURPOSE:
%   Recomputes and updates geometric and pixel statistics for corrected
%   ellipses after intersection refinement. The function:
%     • Updates ellipse area using analytical formula
%     • Rebuilds pixel index lists from ellipse parameters
%     • Falls back to original pixel lists if reconstruction fails
%
% INPUTS:
%   Stats_UNCORRECTED   -> original regionprops structure
%   Stats_CORRECTED     -> corrected regionprops structure (to update)
%   Clean_Image_UpdateS -> reference image (used for sizing)
%   columnsInImage      -> meshgrid columns (not directly used here)
%   rowsInImage         -> meshgrid rows (not directly used here)
%
% OUTPUT:
%   Stats_CORRECTED     -> updated structure with refreshed pixel data
%--------------------------------------------------------------------------

% Final stats generator function
function [Stats_CORRECTED] = FINAL_STATS_FORMULATOR_SMART_INTERSECT_CORRECTOR(varargin)

%--------------------------------------------------------------------------
% Unpack inputs
%--------------------------------------------------------------------------
Stats_UNCORRECTED   = varargin{1};
Stats_CORRECTED     = varargin{2};
Clean_Image_UpdateS = varargin{3};
columnsInImage      = varargin{4}; %#ok<NASGU>
rowsInImage         = varargin{5}; %#ok<NASGU>

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 1: Recompute ellipse areas analytically (πab)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Corrected_Area = (pi * 0.25 .* ...
                 [Stats_CORRECTED(:).MajorAxisLength] .* ...
                 [Stats_CORRECTED(:).MinorAxisLength])';

% Store updated area in structure
for ii = 1:length(Stats_CORRECTED)
    Stats_CORRECTED(ii).Final_Area = Corrected_Area(ii);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 2: Extract ellipse parameters for pixel reconstruction
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

xy_Centroid = vertcat(Stats_CORRECTED.Centroid);

xCenter   = num2cell(xy_Centroid(:,1));
yCenter   = num2cell(xy_Centroid(:,2));
majRadius = num2cell(0.5 * vertcat(Stats_CORRECTED.MajorAxisLength));
minRadius = num2cell(0.5 * vertcat(Stats_CORRECTED.MinorAxisLength));

% Orientation converted to radians (negative to match image coordinates)
phie = num2cell(-deg2rad(vertcat(Stats_CORRECTED.In_plane_Orientation)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 3: Build local bounding grids around each ellipse
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

size_length_2D = size(Clean_Image_UpdateS);

[c, r] = cellfun(@(x,y,ma) ...
    meshgrid(round(x-ma):round(x+ma), ...
             round(y-ma):round(y+ma)), ...
             xCenter, yCenter, majRadius, 'uni', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 4: Determine pixels inside each ellipse (implicit equation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pixel_Idx = cellfun(@(x,y,z,m,n,cg,rg) ...
    find(((cg - x)*cos(z) + (rg - y)*sin(z)).^2/(m)^2 + ...
         (-(cg - x)*sin(z) + (rg - y)*cos(z)).^2/(n)^2 <= 1), ...
         xCenter, yCenter, phie, majRadius, minRadius, c, r, 'uni', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 5: Convert to (x,y) pixel coordinate lists
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pixel_List = cellfun(@(idx,rg,cg) [rg(idx) cg(idx)], ...
                     Pixel_Idx, r, c, 'uni', 0);  % [row col]

% Convert to [x y] ordering
Pixel_List = cellfun(@(x) fliplr(x), Pixel_List, 'uni', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 6: Restore original pixels where reconstruction failed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pixel_List_OLD = {Stats_UNCORRECTED.Pixel_List}';
empty_idx = cell2mat(cellfun(@isempty, Pixel_List, 'uni', 0));

[Pixel_List(empty_idx)] = deal(Pixel_List_OLD(empty_idx));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 7: Constrain pixels to image bounds
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pixel_List = cellfun(@(x) ...
    x((ismember(x(:,1),1:size_length_2D(2)) & ...
       ismember(x(:,2),1:size_length_2D(1))), :), ...
       Pixel_List, 'uni', 0);

% Convert coordinate lists to linear indices
Pixel_Idx = cellfun(@(x) ...
    sub2ind(size_length_2D, x(:,2), x(:,1)), ...
    Pixel_List, 'uni', 0);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP 8: Replace pixel fields in corrected structure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Stats_CORRECTED = rmfield(Stats_CORRECTED, 'Pixel_IDX_List');
Stats_CORRECTED = rmfield(Stats_CORRECTED, 'Pixel_List');

for ii = 1:length(Stats_CORRECTED)
    Stats_CORRECTED(ii).Pixel_IDX_List = Pixel_Idx{ii};
    Stats_CORRECTED(ii).Pixel_List     = Pixel_List{ii};
end