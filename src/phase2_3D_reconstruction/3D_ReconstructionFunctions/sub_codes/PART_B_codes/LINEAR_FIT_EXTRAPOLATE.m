% This function estimates the tortuosity of fiber strands and performs
% directional splitting based on gradient information of fiber sections,
% using a sectioning procedure analogous to a moving centroid approach.
%
% Function workflow:
% (a) Identifies and removes outlier pixel locations using Cook's distance
% (b) Extracts linearly independent points to improve SVD robustness,
%     especially for step-like slopes
% (c) Performs SVD on the filtered points; if a steep/degenerate slope is
%     detected, the original points are preserved for extrapolation

% function [r_3D]      =    LINEAR_FIT_EXTRAPOLATE(varargin)

function [r_3D_cell]      =    LINEAR_FIT_EXTRAPOLATE(varargin)

Pixel_units_fin      =      varargin{1};
blocks_temp          =      varargin{2};
Pixel_units_fin      =      [Pixel_units_fin      (1:size(Pixel_units_fin,1))'];

% extrap_direction      =    varargin{3};
% If extrap_direction   ==  1 → projection upward
% If extrap_direction   ==  2 → projection downward

%% ====>>>> Extract pixels that pass the Cook's distance test
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine rogue (gap) region between bounding fiber blocks
rogue_num                =     blocks_temp{1}(end)+1: blocks_temp{2}(1)-1;
rogue_centroid_region    =     Pixel_units_fin(rogue_num,:);
rogue_depth              =     size(rogue_centroid_region,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP A: If the rogue region is only one slice thick,
% perform simple endpoint extrapolation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(rogue_centroid_region,1) == 1

  r_3D_cell{1}            =      Pixel_units_fin(end,1:2); 
  r_3D_cell{2}            =      Pixel_units_fin(1,1:2);   
  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEP B: Rogue region spans multiple slices → perform linear modeling
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else 
            % Estimate dominant direction of the rogue segment via SVD
            r0                       =     mean(rogue_centroid_region);
            xyz_cent                 =     bsxfun(@minus,rogue_centroid_region,r0);
            [~,~,V]                  =     svd(xyz_cent,0);
            direction_r              =     V(:,1);
            direction_r              =     direction_r ./ vecnorm(direction_r);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            r_3D_cell = cell(1,2);
            
            % Process both extrapolation directions independently
            for ii = 1:2

            extrap_direction    =   ii;
            Centroid_MID_cell_lump = Pixel_units_fin(blocks_temp{ii}(1):blocks_temp{ii}(end),:);

            % If sufficient neighboring samples exist, prefer bounding region fit
            if size(Centroid_MID_cell_lump,1)>= 5  && size(Centroid_MID_cell_lump,1)> size(rogue_centroid_region,1)

                    % Remove outliers using Cook's distance–based regression
                    [centroid,potential_outlier]          =      regstats_func(Centroid_MID_cell_lump);

                    %% Perform SVD on cleaned centroids to obtain fiber direction
                    if ~isempty(centroid)

                            % Keep only linearly independent rows to stabilize SVD
                            tol = 1e-10;
                            [~,rows_idx]                  =      LIROWS(centroid',tol);

                            % If sufficient independent information exists
                            if numel(rows_idx) > 1

                                             r0           =  mean(centroid);
                                             xyz_cent     =  bsxfun(@minus,centroid,r0);
                                             [~,~,V]      =  svd(xyz_cent,0);
                                             direction    =  V(:,1);
                                             direction    =  direction/ vecnorm(direction);

                                    % Guard against degenerate directions
                                    if any(direction == Inf | direction == -Inf | isnan(direction) | direction == -NaN)

                                             % Fallback: straight replication from boundary
                                             if extrap_direction   ==  1
                                             r_3D = repmat(Pixel_units_fin(end,1:2), [rogue_depth,1]); 
                                             else
                                             r_3D = repmat(Pixel_units_fin(1,1:2), [rogue_depth,1]);
                                             end

                                    else                   
                                             % Project points along estimated fiber direction
                                             midpoints = mean(Centroid_MID_cell_lump);
                                             dirvec    = direction';

                                             dl = (rogue_num - midpoints(end))/dirvec(end);                                            
                                             dl = num2cell(dl,1);
                                             points_slice = cellfun(@(x)([midpoints(1) dirvec(1); midpoints(2) dirvec(2)]*[1;x])',dl,'uni',0);
                                             points_slice = vertcat(points_slice{:});
                                             r_3D = points_slice;                
                                    end           

                            else 
                                             % Insufficient rank → fallback to boundary replication
                                             if extrap_direction   ==  1
                                             r_3D = repmat(Pixel_units_fin(end,1:2), [rogue_depth,1]); 
                                             else
                                             r_3D = repmat(Pixel_units_fin(1,1:2), [rogue_depth,1]);
                                             end
                            end

                    else   
                                             % Empty centroid (e.g., uniform Cook's distance)
                                             if extrap_direction   ==  1
                                             r_3D = repmat(Pixel_units_fin(end,1:2), [rogue_depth,1]); 
                                             else
                                             r_3D = repmat(Pixel_units_fin(1,1:2), [rogue_depth,1]);
                                             end
                    end

            % Bounding region exists but is not larger than rogue region
            elseif size(Centroid_MID_cell_lump,1)>= 5  && size(Centroid_MID_cell_lump,1)<= size(rogue_centroid_region,1)
                
                                   if any(direction_r == Inf | direction_r == -Inf | isnan(direction_r) | direction_r == -NaN)

                                             if extrap_direction   ==  1
                                             r_3D = repmat(Pixel_units_fin(end,1:2), [rogue_depth,1]); 
                                             else
                                             r_3D = repmat(Pixel_units_fin(1,1:2), [rogue_depth,1]);
                                             end

                                    else                   
                                             midpoints = mean(Centroid_MID_cell_lump);
                                             dirvec    = direction_r';

                                             dl = (rogue_num - midpoints(end))/dirvec(end);
                                             dl = num2cell(dl,1);
                                             points_slice = cellfun(@(x)([midpoints(1) dirvec(1); midpoints(2) dirvec(2)]*[1;x])',dl,'uni',0);
                                             points_slice = vertcat(points_slice{:});
                                             r_3D = points_slice;                
                                    end           
                
            % Few samples but still larger than rogue region
            elseif  size(Centroid_MID_cell_lump,1)< 5  && size(Centroid_MID_cell_lump,1)> size(rogue_centroid_region,1)
                
                                             r0           =  mean(Centroid_MID_cell_lump);
                                             xyz_cent     =  bsxfun(@minus,Centroid_MID_cell_lump,r0);
                                             [~,~,V]      =  svd(xyz_cent,0);
                                             direction    =  V(:,1);
                                             direction    =  direction ./ vecnorm(direction);

                                    if any(direction == Inf | direction == -Inf | isnan(direction) | direction == -NaN)

                                             if extrap_direction   ==  1
                                             r_3D = repmat(Pixel_units_fin(end,1:2), [rogue_depth,1]); 
                                             else
                                             r_3D = repmat(Pixel_units_fin(1,1:2), [rogue_depth,1]);
                                             end

                                    else                   
                                             midpoints = mean(Centroid_MID_cell_lump);
                                             dirvec    = direction';

                                             dl = (rogue_num - midpoints(end))/dirvec(end);
                                             dl = num2cell(dl,1);
                                             points_slice = cellfun(@(x)([midpoints(1) dirvec(1); midpoints(2) dirvec(2)]*[1;x])',dl,'uni',0);
                                             points_slice = vertcat(points_slice{:});
                                             r_3D = points_slice;                
                                    end           
                
            % Minimal data everywhere → rely on rogue direction
            elseif  size(Centroid_MID_cell_lump,1)< 5  && size(Centroid_MID_cell_lump,1)<= size(rogue_centroid_region,1)
                
                                    if any(direction_r == Inf | direction_r == -Inf | isnan(direction_r) | direction_r == -NaN)

                                             if extrap_direction   ==  1
                                             r_3D = repmat(Pixel_units_fin(end,1:2), [rogue_depth,1]); 
                                             else
                                             r_3D = repmat(Pixel_units_fin(1,1:2), [rogue_depth,1]);
                                             end

                                    else                   
                                             midpoints = mean(Centroid_MID_cell_lump,1);
                                             dirvec    = direction_r';

                                             dl = (rogue_num - midpoints(end))/dirvec(end);
                                             dl = num2cell(dl,1);
                                             points_slice = cellfun(@(x)([midpoints(1) dirvec(1); midpoints(2) dirvec(2)]*[1;x])',dl,'uni',0);
                                             points_slice = vertcat(points_slice{:});
                                             r_3D = points_slice;   
                                    end 
            end
            
            % Store extrapolated segment for this direction
            r_3D_cell{ii} = r_3D;

            end
end

end