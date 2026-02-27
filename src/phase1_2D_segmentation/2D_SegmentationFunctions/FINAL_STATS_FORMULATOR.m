% FINAL_STATS_FORMULATOR
% -------------------------------------------------------------------------
% PURPOSE:
% Recomputes geometric statistics of segmented regions by replacing the
% raw pixel-count area with the area of the fitted bounding ellipse.
% It also reconstructs pixel index lists consistent with the ellipse fit.
%
% INPUTS (via varargin):
%   1) Stats            - Structure containing PixelIdxList for regions
%   2) size_length_2D   - Size of the 2D image slice [rows cols]
%   3) columnsInImage   - Meshgrid of column coordinates (kept for API compatibility)
%   4) rowsInImage      - Meshgrid of row coordinates (kept for API compatibility)
%
% OUTPUT:
%   Stats_FINAL         - Updated statistics structure containing:
%                         • Centroid
%                         • Major/Minor axis lengths
%                         • In-plane orientation
%                         • Out-of-plane orientation
%                         • Final_Area (ellipse-based)
%                         • Pixel_IDX_List
%                         • Pixel_List
% -------------------------------------------------------------------------

function [Stats_FINAL] = FINAL_STATS_FORMULATOR(varargin)

Stats              =  varargin{1};
size_length_2D     =  varargin{2};
columnsInImage     =  varargin{3}; %#ok<NASGU> % retained for compatibility
rowsInImage        =  varargin{4}; %#ok<NASGU> % retained for compatibility

%% ------------------------------------------------------------------------
% Build binary mask from the supplied PixelIdxList
% -------------------------------------------------------------------------
mat_check                                  =     zeros(size_length_2D);
mat_check(vertcat(Stats.PixelIdxList))     =     1;

% Extract ellipse properties from the mask
Stats_FINAL                                =     regionprops(logical(mat_check), ...
                                                'Centroid','MajorAxisLength', ...
                                                'MinorAxisLength','Orientation');

% Rename Orientation → In_plane_Orientation for consistency
[Stats_FINAL(:).In_plane_Orientation]      =     deal(Stats_FINAL(:).Orientation);
Stats_FINAL                                =     rmfield(Stats_FINAL,'Orientation');

%% ------------------------------------------------------------------------
% Compute ellipse-based geometric corrections
% -------------------------------------------------------------------------
% Elliptical area estimate
Corrected_Area = (pi*0.25.* ...
                 [Stats_FINAL(:).MajorAxisLength] .* ...
                 [Stats_FINAL(:).MinorAxisLength])';

% Out-of-plane orientation estimated from axis ratio
Out_of_plane_Orientation = acosd( ...
    [Stats_FINAL(:).MinorAxisLength] ./ ...
    [Stats_FINAL(:).MajorAxisLength]);

% Store computed values
for ii = 1:length(Stats_FINAL)
    Stats_FINAL(ii).Out_of_plane_Orientation = Out_of_plane_Orientation(ii);
    Stats_FINAL(ii).Final_Area               = Corrected_Area(ii);
end

%% ------------------------------------------------------------------------
% Prepare ellipse parameters for pixel reconstruction
% -------------------------------------------------------------------------
xy_Centroid        =    vertcat(Stats_FINAL.Centroid);
xCenter            =    num2cell(xy_Centroid(:,1));
yCenter            =    num2cell(xy_Centroid(:,2));
majRadius          =    num2cell(0.5*vertcat(Stats_FINAL.MajorAxisLength));
minRadius          =    num2cell(0.5*vertcat(Stats_FINAL.MinorAxisLength));

% Convert in-plane orientation to radians (negative sign preserves
% original coordinate convention)
phie               =    num2cell(-deg2rad( ...
                                 vertcat(Stats_FINAL.In_plane_Orientation)));

%% ------------------------------------------------------------------------
% Compute local coordinate grids around each ellipse centroid
% This limits computations to a bounding box per object.
% -------------------------------------------------------------------------
[c,r]= cellfun(@(x,y,ma) ...
        meshgrid(round(x-ma):round(x+ma), ...
                 round(y-ma):round(y+ma)), ...
        xCenter,yCenter,majRadius,'uni',0);

%% ------------------------------------------------------------------------
% Identify pixels that fall inside each rotated ellipse
% -------------------------------------------------------------------------
Pixel_Idx = cellfun(@(x,y,z,m,n,c,r) ...
    find( ((c - x)*cos(z) + (r - y)*sin(z)).^2 /(m)^2  ...
        +(-(c - x)*sin(z) + (r - y)*cos(z)).^2 /(n)^2  <= 1), ...
        xCenter, yCenter, phie, majRadius, minRadius, c, r,'uni',0);

%% ------------------------------------------------------------------------
% Convert local indices to actual pixel coordinate lists
% -------------------------------------------------------------------------
Pixel_List = cellfun(@(x,r,c) [r(x) c(x)], ...
                     Pixel_Idx, r, c,'uni',0);  % [row col]

% Flip to [x y] ordering for consistency with earlier pipeline
Pixel_List = cellfun(@(x)fliplr(x),Pixel_List,'uni',0);

%% ------------------------------------------------------------------------
% Remove coordinates that fall outside image boundaries
% (can occur due to meshgrid expansion near borders)
% -------------------------------------------------------------------------
Pixel_List = cellfun(@(x) ...
    x((ismember(x(:,1),1:size_length_2D(2)) & ...
       ismember(x(:,2),1:size_length_2D(1))),:), ...
    Pixel_List,'uni',0);

%% ------------------------------------------------------------------------
% Convert coordinate lists into linear indices and store results
% -------------------------------------------------------------------------
Pixel_Idx = cellfun(@(x) ...
                    sub2ind(size_length_2D,x(:,2),x(:,1)), ...
                    Pixel_List,'uni',0);

for ii = 1:length(Stats_FINAL)
    Stats_FINAL(ii).Pixel_IDX_List = Pixel_Idx{ii};
    Stats_FINAL(ii).Pixel_List     = Pixel_List{ii};
end

end