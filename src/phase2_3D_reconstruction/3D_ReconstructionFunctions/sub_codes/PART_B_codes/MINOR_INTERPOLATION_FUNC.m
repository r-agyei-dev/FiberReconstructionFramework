function [varargout] = MINOR_INTERPOLATION_FUNC(orig_cent_ext,query_path)

% MINOR_INTERPOLATION_FUNC
% -------------------------------------------------------------------------
% Performs linear interpolation of centroid XY positions across all slice
% indices spanned by the input data.
%
% Purpose:
%   Generates a continuous centroid trajectory by filling in any missing
%   slice locations between the minimum and maximum depth of the provided
%   centroid set.
%
% Method:
%   • Builds a full integer slice vector between min and max depth
%   • Interpolates X and Y coordinates independently using linear interp
%   • Removes one boundary slice depending on scan direction to avoid
%     duplicate overlap with neighboring segments
%
% Inputs:
%   orig_cent_ext - Centroid set [X Y Z] already extrapolated or fitted
%   query_path    - Scan direction flag:
%                   1 → downward scan (remove first slice)
%                   2 → upward scan   (remove last slice)
%
% Output:
%   varargout{1}  - Interpolated centroid trajectory [X Y Z]
%
% Notes:
%   • No extrapolation is performed here (interpolation only)
%   • Assumes slice indices in column 3 are monotonic
%   • Boundary trimming prevents double counting between segments

                        % Build full slice index range to be interpolated
                        query_vec_points = [min(orig_cent_ext(:,3)):max(orig_cent_ext(:,3))];

                        % Known sample slice locations and XY values
                        sample_points_cr = orig_cent_ext(:,3);
                        sample_values_cr = orig_cent_ext(:,1:2);

                        % Interpolate X coordinate along slice direction
                        result(:,1) = interp1(sample_points_cr, ...
                                               sample_values_cr(:,1), ...
                                               query_vec_points', ...
                                               'linear');

                        % Interpolate Y coordinate along slice direction
                        result(:,2) = interp1(sample_points_cr, ...
                                               sample_values_cr(:,2), ...
                                               query_vec_points', ...
                                               'linear');

                        % Reassemble interpolated centroid set
                        orig_cent_extrapolant = [result query_vec_points'];

                        % Remove boundary slice based on scan direction
                        if query_path == 2          % upward scan → drop last slice
                            orig_cent_extrapolant(end,:) = [];
                        elseif query_path == 1      % downward scan → drop first slice
                            orig_cent_extrapolant(1,:) = [];
                        end

                        % Return through varargout for pipeline compatibility
                        varargout{1} = orig_cent_extrapolant;