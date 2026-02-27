function [varargout] = circle_edge_plotter(Data)
%--------------------------------------------------------------------------
% circle_edge_plotter
%
% PURPOSE:
%   Displays the complement of an input image with an overlaid circular
%   boundary and returns a logical mask of the circle pixels.
%
% INPUT:
%   Data - structure containing:
%       .xMeshgrid  -> image width (number of columns)
%       .yMeshgrid  -> image height (number of rows)
%       .image_test -> input grayscale/binary image
%       .radius     -> radius of the circle
%       .xcenter    -> x-coordinate of circle center
%       .ycenter    -> y-coordinate of circle center
%
% OUTPUT:
%   varargout{1} -> logical mask where circle pixels are true
%--------------------------------------------------------------------------

% Extract required fields from input structure
xMeshgrid      =  Data.xMeshgrid;     % Number of columns in the image
yMeshgrid      =  Data.yMeshgrid;     % Number of rows in the image
image_test     =  Data.image_test;    % Image to display
radius         =  Data.radius;        % Circle radius
xcenter        =  Data.xcenter;       % Circle center (x)
ycenter        =  Data.ycenter;       % Circle center (y)

%----------------------------------------------------------------------
% Create coordinate grids for pixel indexing
% columnsInImage -> x-coordinates
% rowsInImage    -> y-coordinates
%----------------------------------------------------------------------
[columnsInImage,rowsInImage] = meshgrid(1:xMeshgrid, 1:yMeshgrid);

%----------------------------------------------------------------------
% Display the complemented image and overlay the circle boundary
%----------------------------------------------------------------------
figure,imshow(imcomplement(image_test))  % Show inverted image
hold on  

% Compute top-left corner of bounding box for the circle
xs = xcenter - radius; 
ys = ycenter - radius;

% Rectangle position vector: [x y width height]
pos = [xs ys 2*radius 2*radius];

% Draw circular boundary using rectangle with full curvature
rectangle('Position',pos,'Curvature',[1 1], ...
          'LineWidth',3,'EdgeColor','g');

hold off 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------------------------------------------------------
% Create logical mask for pixels inside the circle
% Equation of circle: (x - xc)^2 + (y - yc)^2 <= r^2
%----------------------------------------------------------------------
circlePixels = (rowsInImage - ycenter).^2 + ...
               (columnsInImage - xcenter).^2 <= radius.^2;

%----------------------------------------------------------------------
% Return output
%----------------------------------------------------------------------
varargout{1} = circlePixels;

end