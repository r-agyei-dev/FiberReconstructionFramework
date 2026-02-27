function [xyq, grad] = CENTROID_STAGGER_CORRRECTOR_DIRECT(Centroid_MID_cell_lump, fig_switch, size_length)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% CENTROID_STAGGER_CORRRECTOR_DIRECT
%
% This function refits the centroids of a fiber in 3D space to correct for any
% staggered misalignment along the stack. It performs a 3D linear fit of the 
% centroid points and constrains points to stay within the bounds of the volume.
%
% Inputs:
% - Centroid_MID_cell_lump : Nx2 or Nx3 array of centroids for a single fiber segment
% - fig_switch             : 0 (no plots) or 1 (plot the original vs. refitted centroids)
% - size_length            : 1x2 or 1x3 array specifying the maximum allowable dimensions
%
% Outputs:
% - xyq  : Refitted centroid positions constrained within the volume
% - grad : Gradient of the fiber along the 3D fit (approximate orientation)
%
% NOTE:
% - The number of centroids must be greater than 2; otherwise MATLAB regstats 
%   or linear fitting will fail with "more predictor variables than observations."
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if size(Centroid_MID_cell_lump, 1) > 2
    % Perform a 3D linear fit through the centroids to correct staggered alignment
    Pixel_units_fin = Centroid_MID_cell_lump;
    [r_3D, grad] = LINEAR_FIT_FUNCTION_3D_WITH_FIG_DIRECT(Pixel_units_fin, fig_switch);

    % Constrain centroids that exceed the volume bounds to the edge values
    r_3D(r_3D(:,1) > size_length(2), 1) = size_length(2);
    r_3D(r_3D(:,2) > size_length(1), 2) = size_length(1);

    xyq = r_3D;

    % Optional 3D scatter plot for debugging
    if fig_switch == 1
        figure, scatter3(Pixel_units_fin(:,1), Pixel_units_fin(:,2), (1:size(Pixel_units_fin,1))')
        hold on
        scatter3(xyq(:,1), xyq(:,2), (1:size(xyq,1))','r'), grid on
    end

else
    % If only 1 or 2 points, perform a simplified correction to keep sizes consistent
    if size(Centroid_MID_cell_lump,2) == 2
        xyq = [Centroid_MID_cell_lump, (1:size(Centroid_MID_cell_lump,1))'];
        % Approximate gradient based on the line from first to last centroid
        grad = min([180 - acosd(dot(normr(xyq(1,:) - xyq(end,:)), [0 0 1])); ...
                     acosd(dot(normr(xyq(1,:) - xyq(end,:)), [0 0 1]))]);
    else
        xyq = Centroid_MID_cell_lump;
        grad = 0;
    end
end