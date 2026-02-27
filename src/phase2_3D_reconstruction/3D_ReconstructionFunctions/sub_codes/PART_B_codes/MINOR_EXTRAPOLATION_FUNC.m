% MINOR_EXTRAPOLATION_FUNC
% -------------------------------------------------------------------------
% Performs linear interpolation/extrapolation of centroid XY coordinates
% along the slice (Z) direction.
%
% Purpose:
%   Given a fitted centroid trajectory and the original centroid samples,
%   this function reconstructs (or extrapolates) the XY centroid positions
%   at the slice locations defined by the fitted data.
%
% Method:
%   • Uses 1D linear interpolation along the third column (slice index)
%   • Applies extrapolation outside the sampled range when needed
%   • Processes X and Y coordinates independently
%
% Inputs:
%   fitted_cent - Target centroid positions with desired slice indices
%                 (must contain column 3 as the query depth)
%   orig_cent   - Original centroid samples used as interpolation source
%                 (columns 1–2 = XY, column 3 = slice index)
%
% Output:
%   varargout{1} - Extrapolated/interpolated centroid set in the format
%                  [X Y Z], aligned to fitted_cent slice locations
%
% Notes:
%   • Linear extrapolation is enabled via 'extrap'
%   • Assumes slice index (column 3) is monotonic in orig_cent

function [varargout] = MINOR_EXTRAPOLATION_FUNC(fitted_cent,orig_cent)

                        % Target slice locations where centroids are needed
                        query_vec_points = fitted_cent(:,3);

                        % Known sample slice locations and XY values
                        sample_points_cr = orig_cent(:,3);
                        sample_values_cr = orig_cent(:,1:2);

                        % Interpolate/extrapolate X coordinate
                        result(:,1) = interp1(sample_points_cr', ...
                                               sample_values_cr(:,1), ...
                                               query_vec_points', ...
                                               'linear','extrap');

                        % Interpolate/extrapolate Y coordinate
                        result(:,2) = interp1(sample_points_cr', ...
                                               sample_values_cr(:,2), ...
                                               query_vec_points', ...
                                               'linear','extrap');

                        % Reassemble full centroid with slice index
                        orig_cent_extrapolant = [result fitted_cent(:,3)];

                        % Return through varargout for pipeline compatibility
                        varargout{1} = orig_cent_extrapolant;