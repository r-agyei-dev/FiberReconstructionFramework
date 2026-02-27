function [Stats_FINAL]=FINAL_STATS_FORMULATOR_SMART_SH_TT(varargin)
% FINAL_STATS_FORMULATOR_SMART_SH_TT
% -------------------------------------------------------------------------
% PURPOSE:
% Builds a refined statistics structure from connected-component data.
% It computes ellipse-based geometric properties and reconstructs the
% pixel lists corresponding to each detected object.
%
% INPUTS (via varargin):
%   1) Stats              - Structure containing PixelIdxList for objects
%   2) size_length_2D     - Size of the 2D image [rows cols]
%   3) columnsInImage     - Full column meshgrid (currently unused)
%   4) rowsInImage        - Full row meshgrid (currently unused)
%
% OUTPUT:
%   Stats_FINAL           - Updated statistics structure with:
%                           • Centroid
%                           • Major/Minor axis lengths
%                           • In-plane orientation
%                           • Out-of-plane orientation
%                           • Final_Area (ellipse-based)
%                           • Pixel_IDX_List
%                           • Pixel_List
% -------------------------------------------------------------------------

Stats              =  varargin{1};
size_length_2D     =  varargin{2};
columnsInImage     =  varargin{3}; %#ok<NASGU> % retained for compatibility
rowsInImage        =  varargin{4}; %#ok<NASGU> % retained for compatibility

%% ------------------------------------------------------------------------
% Build a binary mask from the provided PixelIdxList
% -------------------------------------------------------------------------
mat_check                                  =     zeros(size_length_2D);
mat_check(vertcat(Stats.PixelIdxList))     =     1;

% Extract ellipse statistics from the binary mask
Stats_FINAL                                =     regionprops(logical(mat_check), ...
                                                'Centroid','MajorAxisLength', ...
                                                'MinorAxisLength','Orientation');

% Rename Orientation to In_plane_Orientation for clarity
[Stats_FINAL(:).In_plane_Orientation]      =     deal(Stats_FINAL(:).Orientation);
Stats_FINAL                                =     rmfield(Stats_FINAL,'Orientation');

%% ------------------------------------------------------------------------
% Compute corrected geometric quantities
% -------------------------------------------------------------------------
% Elliptical area approximation
Corrected_Area = (pi*0.25.* ...
                 [Stats_FINAL(:).MajorAxisLength] .* ...
                 [Stats_FINAL(:).MinorAxisLength])';

% Out-of-plane tilt estimate from axis ratio
Out_of_plane_Orientation = acosd( ...
    [Stats_FINAL(:).MinorAxisLength] ./ ...
    [Stats_FINAL(:).MajorAxisLength]);

% Store computed quantities
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
% previous coordinate convention)
phie               =    num2cell(-deg2rad( ...
                                 vertcat(Stats_FINAL.In_plane_Orientation)));

%% ------------------------------------------------------------------------
% Generate local meshgrids around each ellipse
% This limits computations to a bounding box per object (efficient)
% -------------------------------------------------------------------------
[c,r]= cellfun(@(x,y,ma) ...
        meshgrid(round(x-ma):round(x+ma), ...
                 round(y-ma):round(y+ma)), ...
        xCenter,yCenter,majRadius,'uni',0);

%% ------------------------------------------------------------------------
% Determine which pixels fall inside each rotated ellipse
% (Pseudo-linear index inside the local grid)
% -------------------------------------------------------------------------
Pixel_Idx = cellfun(@(x,y,z,m,n,c,r) ...
    find( ((c - x)*cos(z) + (r - y)*sin(z)).^2 /(m)^2  ...
        +(-(c - x)*sin(z) + (r - y)*cos(z)).^2 /(n)^2  <= 1), ...
        xCenter, yCenter, phie, majRadius, minRadius, c, r, 'uni', 0);

%% ------------------------------------------------------------------------
% Convert local indices into actual (x,y) pixel coordinate lists
% -------------------------------------------------------------------------
Pixel_List = cellfun(@(x,r,c) [r(x) c(x)], ...
                     Pixel_Idx, r, c, 'uni', 0);   % [row col]

% Flip to [x y] ordering for consistency with earlier code
Pixel_List = cellfun(@(x)fliplr(x), Pixel_List, 'uni', 0);

%% ------------------------------------------------------------------------
% Constrain pixel coordinates so they stay inside image bounds
% -------------------------------------------------------------------------
Pixel_List = cellfun(@(x) ...
    x((ismember(x(:,1),1:size_length_2D(2)) & ...
       ismember(x(:,2),1:size_length_2D(1))), :), ...
    Pixel_List,'uni',0);

%% ------------------------------------------------------------------------
% Convert coordinate lists to linear indices
% -------------------------------------------------------------------------
Pixel_Idx = cellfun(@(x) ...
                    sub2ind(size_length_2D,x(:,2),x(:,1)), ...
                    Pixel_List,'uni',0);

%% ------------------------------------------------------------------------
% Store results back into Stats_FINAL
% -------------------------------------------------------------------------
for ii = 1:length(Stats_FINAL)
    Stats_FINAL(ii).Pixel_IDX_List = Pixel_Idx{ii};
    Stats_FINAL(ii).Pixel_List     = Pixel_List{ii};
end

end