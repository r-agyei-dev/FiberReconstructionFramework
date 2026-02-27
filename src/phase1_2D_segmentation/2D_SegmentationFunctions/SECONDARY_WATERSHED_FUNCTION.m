%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SECONDARY_WATERSHED_FUNCTION
% -------------------------------------------------------------------------
% PURPOSE:
% This function performs a refined watershed segmentation on a 2D image, 
% computes geometrical metadata of the resulting ellipses, and calculates 
% additional statistics such as in-plane and out-of-plane orientation.
%
% INPUTS:
%   Image_crop - Grayscale 2D image slice to be segmented
%   level      - Threshold for binarization (Otsu or manual)
%   mask_value - Parameter controlling minima suppression for watershed
%
% OUTPUTS:
%   Stats - Structure array containing segmented region metadata:
%           Centroid, MajorAxisLength, MinorAxisLength, In_plane_Orientation,
%           Out_of_plane_Orientation, Pixel_IDX_List, Pixel_List, Final_Area
%
% DESCRIPTION:
%   1. Converts the image to binary using the specified threshold.
%   2. Computes the distance transform for watershed segmentation.
%   3. Applies minima suppression with mask_value and imposes it on the distance transform.
%   4. Watershed segmentation separates touching regions.
%   5. Computes ellipse-based metrics for each region (area, orientation, pixel coordinates).
% -------------------------------------------------------------------------

function [Stats] = SECONDARY_WATERSHED_FUNCTION(varargin)

%% ------------------------------------------------------------------------
% Step 1: Unpack input arguments
% -------------------------------------------------------------------------
Image_crop = varargin{1};
level      = varargin{2};
mask_value = varargin{3};

%% ------------------------------------------------------------------------
% Step 2: Binarize image using threshold
% -------------------------------------------------------------------------
Image_crop = imbinarize(Image_crop, level);  % Convert grayscale to binary

%% ------------------------------------------------------------------------
% Step 3: Distance transform and watershed
% -------------------------------------------------------------------------
D  = -bwdist(~Image_crop, 'euclidean');  % Negative distance transform
Ld2 = watershed(imimposemin(D, imextendedmin(D, mask_value)));  % Impose minima and segment
Image_crop(Ld2 == 0) = 0;  % Set watershed ridge lines to background

%% ------------------------------------------------------------------------
% Step 4: Extract region properties
% -------------------------------------------------------------------------
Stats = regionprops(bwareaopen(Image_crop,5), ...
    'Centroid','MinorAxisLength','MajorAxisLength','PixelIdxList','Orientation'); 
% Removes small objects (<5 pixels) and returns geometrical features

%% ------------------------------------------------------------------------
% Step 5: Compute in-plane and out-of-plane orientation
% -------------------------------------------------------------------------
[Stats(:).In_plane_Orientation] = deal(Stats(:).Orientation);  % Copy orientation
Stats = rmfield(Stats,'Orientation');  % Remove redundant field
Stats = rmfield(Stats,'PixelIdxList'); % Will recompute later as Pixel_List

Corrected_Area           =  (pi*0.25 .* [Stats(:).MajorAxisLength] .* [Stats(:).MinorAxisLength])';  % Area of ellipse
Out_of_plane_Orientation =  acosd([Stats(:).MinorAxisLength] ./ [Stats(:).MajorAxisLength]);          % Angle from minor/major ratio

for ii = 1:length(Stats)
    Stats(ii).Out_of_plane_Orientation = Out_of_plane_Orientation(ii);
    Stats(ii).Final_Area               = Corrected_Area(ii);  % Update area
end

%% ------------------------------------------------------------------------
% Step 6: Compute pixel coordinates for each ellipse
% -------------------------------------------------------------------------
xy_Centroid  = vertcat(Stats.Centroid);
xCenter      = num2cell(xy_Centroid(:,1));
yCenter      = num2cell(xy_Centroid(:,2));
majRadius    = num2cell(0.5*vertcat(Stats.MajorAxisLength));
minRadius    = num2cell(0.5*vertcat(Stats.MinorAxisLength));
phie         = num2cell(-deg2rad(vertcat(Stats.In_plane_Orientation)));

size_length_2D = size(Image_crop);

% Generate meshgrid around centroid for candidate pixels
[c,r] = cellfun(@(x,y,ma) meshgrid(round(x-ma):round(x+ma), round(y-ma):round(y+ma)), ...
    xCenter, yCenter, majRadius, 'uni', 0);

% Compute linear indices of pixels within ellipse boundary
Pixel_Idx = cellfun(@(x,y,z,m,n,c,r) find( ((c - x)*cos(z) + (r - y)*sin(z)).^2 /(m)^2  ...
                                     + (-(c - x)*sin(z) + (r - y)*cos(z)).^2 /(n)^2 <= 1), ...
                     xCenter, yCenter, phie, majRadius, minRadius, c, r, 'uni', 0);

% Convert linear indices to row-column lists
Pixel_List = cellfun(@(x,r,c) [r(x) c(x)], Pixel_Idx, r, c, 'uni', 0); 
Pixel_List = cellfun(@(x) fliplr(x), Pixel_List, 'uni', 0);  % Flip to [col,row]

% Remove any out-of-bound coordinates
Pixel_List = cellfun(@(x) x((ismember(x(:,1), 1:size_length_2D(2)) & ismember(x(:,2), 1:size_length_2D(1))), :), ...
                     Pixel_List, 'uni', 0);

% Add linear indices and pixel lists to metadata
Pixel_Idx = cellfun(@(x) sub2ind(size_length_2D, x(:,2), x(:,1)), Pixel_List, 'uni', 0);
for ii = 1:length(Stats)
    Stats(ii).Pixel_IDX_List = Pixel_Idx{ii};
    Stats(ii).Pixel_List     = Pixel_List{ii};
end

end