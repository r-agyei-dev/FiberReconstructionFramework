function [xe,ye] = ellipse_plot_func_circ_post_process(Stats,matrix_Mid,new_idx,thresh_max,thresh_min)
%--------------------------------------------------------------------------
% ellipse_plot_func_circ
%
% PURPOSE:
%   Overlays fitted ellipses (from regionprops-style metadata) onto an
%   image and returns the ellipse boundary coordinates.
%
% INPUTS:
%   Stats       -> structure array containing ellipse parameters:
%                  .Centroid
%                  .MajorAxisLength
%                  .MinorAxisLength
%                  .In_plane_Orientation (degrees)
%   matrix_Mid  -> image on which ellipses are drawn
%   new_idx     -> index used only for display title
%
% OUTPUTS:
%   xe, ye      -> cell arrays containing x and y ellipse coordinates
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% Display the image that will receive the ellipse overlays
%--------------------------------------------------------------------------
% figure, imshow(matrix_Mid), title(sprintf('Image%d',new_idx))
figure
subplot(1,2,1)          % Create a subplot
imshow(matrix_Mid, [])   % Display the image in that subplot
title(sprintf('Image %d', new_idx))  % Add a title

subplot(1,2,2)          % Create a subplot
imshow(matrix_Mid, [])   % Display the image in that subplot
title(sprintf('Image %d', new_idx))  % Add a title
hold on

%--------------------------------------------------------------------------
% Parameterize a unit ellipse (used as template for all regions)
%--------------------------------------------------------------------------
phie   = linspace(0,2*pi,500);   % angle samples for smooth ellipse
cosphi = cos(phie);
sinphi = sin(phie);



Stats = Stats(vertcat(Stats.MajorAxisLength)>thresh_max);
Stats = Stats(vertcat(Stats.MinorAxisLength)>(thresh_min));
%--------------------------------------------------------------------------
% Loop through each detected region and draw its ellipse
%--------------------------------------------------------------------------
for k = 1:length(Stats)

    % Extract ellipse center
    xbar = Stats(k).Centroid(1);
    ybar = Stats(k).Centroid(2);

    % Semi-axis lengths (regionprops gives full lengths)
    a = Stats(k).MajorAxisLength/2;
    b = Stats(k).MinorAxisLength/2;
    
    % Convert orientation from degrees to radians
    thetae = pi*Stats(k).In_plane_Orientation/180;

    % Rotation matrix for ellipse orientation
    R = [ cos(thetae)   sin(thetae)
         -sin(thetae)   cos(thetae)];
    
    % Generate ellipse in local coordinates
    xy = [a*cosphi; b*sinphi];

    % Rotate ellipse to match region orientation
    xy = R*xy;
    
    % Translate ellipse to its centroid location
    xe{k} = xy(1,:) + xbar;
    ye{k} = xy(2,:) + ybar;
  
    % Plot ellipse boundary
    plot(xe{k},ye{k},'r','LineWidth',1);
end

hold off
end