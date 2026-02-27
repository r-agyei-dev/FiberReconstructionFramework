function [varargout] = ALIGNMENT_FUNCTION(varargin)

% INPUTS:
% total_centroid_temp    -> all centroid coordinates of the fiber pixels
% bounding_ellipse_index -> start and end indices of ellipses along the fiber

% OUTPUT:
% varargout{1}           -> euclidean deviations between projected and original centroids

total_centroid_temp    = varargin{1};
bounding_ellipse_index = varargin{2};

% Initialize array to store deviations
euclidean_deviationpp = zeros(1, size(bounding_ellipse_index, 1));

% Loop over each segment defined by bounding ellipses
for mm = 1:size(bounding_ellipse_index, 1)
    
    % Determine indices for the current segment
    if mm == 1
        var             = 1:bounding_ellipse_index(1,1);
        Pixel_units_fin = total_centroid_temp(var,:);
        extra_slice     = bounding_ellipse_index(1,2) - var(end);
    else 
        var             = bounding_ellipse_index(mm-1,2):bounding_ellipse_index(mm,1); 
        Pixel_units_fin = total_centroid_temp(var,:);
        extra_slice     = bounding_ellipse_index(mm,2) - var(end);
    end
    
    % ======================= FIBER ALIGNMENT & 3D PROJECTION ========================
    % This section computes a 3D projected line for the fiber segment using centroid data:
    % 1. Convert pixel locations into 3D points along the fiber axis.
    % 2. Use singular value decomposition (SVD) to determine the principal direction of the fiber.
    % 3. Correct the fiber line for steep slopes or straight-line cases.
    % 4. Include extra slices if needed for alignment.

    Centroid_MID_cell_lump = Pixel_units_fin;

    if size(Centroid_MID_cell_lump,1) < 2
        % If only one point, construct a trivial 3D centroid with slice index
        Centroid_MID_cell_lump_3D = [Centroid_MID_cell_lump, (1:size(Centroid_MID_cell_lump,1))'];
        r_3D = Centroid_MID_cell_lump_3D;
    else
        % Add slice index to each centroid to create 3D representation
        Centroid_MID_cell_lump_3D = [Centroid_MID_cell_lump, (1:size(Centroid_MID_cell_lump,1))'];   

        % Compute principal direction using SVD
        centroid = Centroid_MID_cell_lump_3D;
        r0       = mean(centroid);
        xyz_cent = bsxfun(@minus, centroid, r0);
        [~,~,V]  = svd(xyz_cent, 0);
        direction = V(:,1);
        direction = normc(direction);

        % ===================== STRAIGHT LINE CORRECTION =========================
        % Handle edge cases where direction contains Inf, NaN, or invalid values
        if any(direction == Inf | direction == -Inf | isnan(direction) | direction == -NaN)
            r_3D = [repmat(centroid(1,1:2), [size(Centroid_MID_cell_lump,1),1]) ...
                     (1:size(Centroid_MID_cell_lump,1)+extra_slice)'];
            centroid_3D = Centroid_MID_cell_lump_3D;
        else
            % Project centroids along principal direction to generate 3D fiber
            midpoints = mean(Centroid_MID_cell_lump_3D);
            dirvec    = direction';

            slice_num = 1:size(Centroid_MID_cell_lump_3D,1)+extra_slice;
            dl        = (slice_num - midpoints(end)) / dirvec(end);
            dl        = num2cell(dl,1);

            points_slice = cellfun(@(x)([midpoints(1) dirvec(1); midpoints(2) dirvec(2)]*[1;x])', dl, 'uni', 0);
            points_slice = vertcat(points_slice{:});
            points_slice(:,end+1) = (1:size(points_slice,1))';
            r_3D = points_slice;                
        end
    end

    % Compute deviation between projected fiber and original last centroid
    euclidean_deviationpp(mm) = norm(r_3D(end,1:2) - total_centroid_temp(bounding_ellipse_index(mm,2),:));
end                            

% Return euclidean deviations
varargout{1} = euclidean_deviationpp;

end