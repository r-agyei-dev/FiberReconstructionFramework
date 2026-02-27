%% ALIGNMENT_FUNCTION2
% This function determines whether "rogue" regions (fiber segments between bounding ellipses)
% should be extrapolated or stitched based on the in-plane deviation of their centroids.

function [varargout] = ALIGNMENT_FUNCTION2(varargin)

% INPUTS:
% total_centroid_temp    -> matrix of centroid coordinates for all pixels
% bounding_ellipse_index -> start and end indices of bounding ellipses
% blocks                 -> optional blocks of fiber segments (not used in this section)

% OUTPUT:
% varargout{1}           -> binary array indicating if rogue regions need extrapolation (1) or not (0)

total_centroid_temp    = varargin{1};
bounding_ellipse_index = varargin{2};
blocks                 = varargin{3};

% Initialize array to store deviation/extrapolation flags
euclidean_deviation = zeros(1, size(bounding_ellipse_index,1));
tol_var = 2.5;  % distance tolerance to decide if rogue region should be processed

% Loop through each pair of bounding ellipses
for tt = 1:size(bounding_ellipse_index,1)

    % Compute distances between bounding ellipse centroids and their immediate neighbors
    new_dist_1 = sqrt(sum((total_centroid_temp(bounding_ellipse_index(tt,1),:) - ...
                           total_centroid_temp(bounding_ellipse_index(tt,1)+1,:)).^2));    
    new_dist_2 = sqrt(sum((total_centroid_temp(bounding_ellipse_index(tt,2),:) - ...
                           total_centroid_temp(bounding_ellipse_index(tt,2)-1,:)).^2));    
  
    % Only consider extrapolation if both distances exceed tolerance
    if ~(new_dist_1 <= tol_var || new_dist_2 <= tol_var)
        
        % Compute the mean centroid of the rogue region between bounding ellipses
        mean_centroid = mean(total_centroid_temp(bounding_ellipse_index(tt,1)+1 : ...
                                                 bounding_ellipse_index(tt,2)-1,:), 1);
        
        % Check if the rogue region is located within the bounding ellipse region
        if (mean_centroid(1) > min(total_centroid_temp(bounding_ellipse_index(tt,1),1), ...
                                    total_centroid_temp(bounding_ellipse_index(tt,2),1))) && ...
           (mean_centroid(1) < max(total_centroid_temp(bounding_ellipse_index(tt,1),1), ...
                                    total_centroid_temp(bounding_ellipse_index(tt,2),1))) && ...   
           (mean_centroid(2) > min(total_centroid_temp(bounding_ellipse_index(tt,1),2), ...
                                    total_centroid_temp(bounding_ellipse_index(tt,2),2))) && ...
           (mean_centroid(2) < max(total_centroid_temp(bounding_ellipse_index(tt,1),2), ...
                                    total_centroid_temp(bounding_ellipse_index(tt,2),2)))
        
            % Flag for extrapolation
            euclidean_deviation(tt) = 1;   
        end   
    end
end

% Return the deviation/extrapolation flag array
varargout{1} = euclidean_deviation;

end