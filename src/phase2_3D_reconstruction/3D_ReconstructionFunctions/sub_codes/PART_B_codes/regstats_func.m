function [fin_out,non_potential_outlier] = regstats_func(Centroid_MID_cell_lump)

% REGSTATS_FUNC
% -------------------------------------------------------------------------
% Performs a simple linear regression on the centroid trajectory and
% removes points that behave like outliers based on Cook's distance.
%
% INPUT
%   Centroid_MID_cell_lump : Nx2 matrix of centroid coordinates [x y]
%
% OUTPUTS
%   fin_out                 : filtered centroid coordinates (outliers removed)
%   non_potential_outlier   : indices of points retained after filtering
% -------------------------------------------------------------------------

% Flip columns so the regression is performed in the expected orientation.
% (Original data assumed [x y]; regstats uses Y vs X.)
xy  =  fliplr(Centroid_MID_cell_lump);

% Optional visualization for debugging (disabled)
% figure, scatter(xy(:,end),xy(:,1),'*')

% Perform linear regression:
%   Y = xy(:,end)
%   X = xy(:,1)
stats = regstats(xy(:,end),xy(:,1),'linear');

% NOTE:
% Cook's distance is used to identify influential points.
% A common rule-of-thumb threshold is n/4, but here we instead
% keep points within ±1 standard deviation of the mean Cook's distance.

% Indices of points considered NON-outliers
non_potential_outlier = find( ...
    (mean(stats.cookd) - std(stats.cookd)) < stats.cookd & ...
     stats.cookd < (mean(stats.cookd) + std(stats.cookd)));

% Return the filtered centroid coordinates in original orientation
fin_out = fliplr(xy(non_potential_outlier,:));

% -------------------------------------------------------------------------
% Additional debugging / experimental code (kept commented intentionally)
% -------------------------------------------------------------------------
% figure, scatter(xy(potential_outlier,end),xy(potential_outlier,1),'*')
% A_mat = [xy ones(length(xy),1)];
% B_mat = (1:size(xy,1))';
% g = A_mat\B_mat