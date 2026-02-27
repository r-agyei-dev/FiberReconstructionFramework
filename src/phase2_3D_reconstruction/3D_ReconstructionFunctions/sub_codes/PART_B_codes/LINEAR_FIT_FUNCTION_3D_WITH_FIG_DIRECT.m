% This function estimates the local tortuosity and principal direction of
% a fiber strand using a moving-centroid style linear model.
%
% Processing steps:
% (a) Removes outlier pixel locations using Cook's distance
% (b) Retains only linearly independent samples to stabilize SVD,
%     especially for step-like or near-vertical slopes
% (c) Uses SVD to estimate the dominant fiber direction
% (d) Generates a fitted 3D centerline and computes the inclination angle

function [r_3D,grad] = LINEAR_FIT_FUNCTION_3D_WITH_FIG_DIRECT(Pixel_units_fin,fig_switch)

% Treat input as the working centroid set
Centroid_MID_cell_lump = Pixel_units_fin;

% Ensure a 3D representation by appending slice index if only XY provided
if size(Centroid_MID_cell_lump,2) == 2
    Centroid_MID_cell_lump_3D = [Centroid_MID_cell_lump, (1:size(Centroid_MID_cell_lump,1))'];
else
    Centroid_MID_cell_lump_3D = Centroid_MID_cell_lump;   
end

% Remove outliers via Cook's distance regression
[centroid,potential_outlier] = regstats_func(Centroid_MID_cell_lump);

if ~isempty(centroid)

    % Keep only linearly independent rows for numerical stability
    tol = 1e-10;
    [~,rows_idx] = LIROWS(centroid',tol);

    % Append outlier flag (or auxiliary info) to centroid set
    centroid(:,end+1) = potential_outlier;
             
    % Proceed only if sufficient independent information exists
    if numel(rows_idx) > 1

        % ====>>>> Estimate dominant direction using SVD
        r0       = mean(centroid);
        xyz_cent = bsxfun(@minus,centroid,r0);
        [~,~,V]  = svd(xyz_cent,0);
        direction = V(:,1);
        direction = normc(direction);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step A: Direct parametric line reconstruction

        % Guard against degenerate or invalid directions
        if any(direction == Inf | direction == -Inf | isnan(direction) | direction == -NaN)

            % Fallback: replicate first centroid across slices
            r_3D        = [repmat(centroid(1,1:2), [size(Centroid_MID_cell_lump,1),1]) ...
                           (1:size(Centroid_MID_cell_lump,1))'];
            centroid_3D = Centroid_MID_cell_lump_3D;

        else                   
            % Project along principal direction
            midpoints = mean(Centroid_MID_cell_lump_3D);
            dirvec    = direction';
                         
            slice_num = 1 : size(Centroid_MID_cell_lump_3D,1);
            dl        = (slice_num - midpoints(end)) / dirvec(end);
                         
            dl           = num2cell(dl,1);
            points_slice = cellfun(@(x)([midpoints(1) dirvec(1); ...
                                         midpoints(2) dirvec(2)]*[1;x])', ...
                                         dl,'uni',0);
            points_slice = vertcat(points_slice{:});
            points_slice(:,end+1) = (1:size(points_slice,1))';
            r_3D = points_slice;                
        end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % ====>>>> Insufficient independent rows → return original data
    else 
        centroid_3D = Centroid_MID_cell_lump_3D;
        r_3D        = Centroid_MID_cell_lump_3D;
    end

% ====>>>> Cook's distance removed everything → fallback to original data
else   
    centroid_3D = Centroid_MID_cell_lump_3D;
    r_3D        = Centroid_MID_cell_lump_3D;
end

%% Optional visualization for verification
if fig_switch == 1
    figure, scatter3(centroid_3D(:,1),centroid_3D(:,2),centroid_3D(:,3),'b','LineWidth',3)
    hold on 
    scatter3(r_3D(:,1),r_3D(:,2),r_3D(:,3),'r','filled'), grid on
    hold off
end

%% Compute inclination angle relative to the Z-axis
direct_vect = normr(r_3D(1,:) - r_3D(end,:));
grad = min([180 - acosd(dot(direct_vect,[0 0 1])); ...
            acosd(dot(direct_vect,[0 0 1]))]);

end