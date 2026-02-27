%% =========================================================================
% PURPOSE:
% This function performs extrapolation for "rogue" (missing) centroid
% regions at the beginning or end of a fiber track.
%
% It estimates the missing centroid positions by fitting a direction
% vector using SVD and projecting the trajectory forward/backward.
%
% INPUTS:
%   non_rogue_area_find → logical vector indicating valid (1) vs rogue (0)
%   Interpolant_Result  → interpolated centroid matrix (Nx2)
%   total_centroid_temp → original centroid coordinates
%   size_length         → image size [rows cols]
%
% OUTPUT:
%   varargout{1} → corrected centroid coordinates (Nx2)
%% =========================================================================

function [varargout] = rogue_extrapolation(varargin)

   % ----------------------------------------------------------------------
   % Unpack inputs
   % ----------------------------------------------------------------------
   non_rogue_area_find    =  varargin{1};
   Interpolant_Result     =  varargin{2};
   total_centroid_temp    =  varargin{3};
   size_length            =  varargin{4};
   
   % ----------------------------------------------------------------------
   % Append slice index as third column to centroid matrix
   % (creates [x y slice#])
   % ----------------------------------------------------------------------
   total_centroid_TMP   = [total_centroid_temp (1:size(total_centroid_temp,1))'];

   % ----------------------------------------------------------------------
   % Check if rogue regions exist at the beginning or end
   % ----------------------------------------------------------------------
   if (ismember(1,find(non_rogue_area_find == 0)) || ...
       ismember(numel(non_rogue_area_find),find(non_rogue_area_find == 0)))
    
        % Locate contiguous rogue blocks
        [blocks] = zeros_function_locater(~non_rogue_area_find);
    
        %% =================================================================
        % CASE 1: Rogue region at the START of the fiber
        % =================================================================
        if ismember(1,find(non_rogue_area_find == 0))
            
            % Extract rogue centroid block
            rogue_volume   = total_centroid_TMP(blocks{1},:);
            rogue_num_orig = rogue_volume(:,3);
            
            % --------------------------------------------------------------
            % Remove outliers if multiple rogue points exist
            % --------------------------------------------------------------
            if numel(rogue_num_orig) > 1
                [rogue_volume,~] = regstats_func(rogue_volume);
                rogue_num_curr   = rogue_volume(:,3);
            else 
                rogue_num_curr   = rogue_num_orig;
            end
            
            % --------------------------------------------------------------
            % Fallback: insufficient data for direction estimation
            % --------------------------------------------------------------
            if (numel(rogue_num_curr) == 1 || isempty(rogue_num_curr))
                
                r_3D = repmat(total_centroid_TMP(rogue_num_orig(end)+1,1:2), ...
                              [numel(rogue_num_orig),1]); 
                Interpolant_Result(blocks{1},1:2) = r_3D;
                
            else 
                % ----------------------------------------------------------
                % Estimate principal direction using SVD
                % ----------------------------------------------------------
                r0        = mean(rogue_volume);
                xyz_cent  = bsxfun(@minus,rogue_volume,r0);
                [~,~,V]   = svd(xyz_cent,0);
                direction_r = V(:,1);
                direction_r = direction_r ./ vecnorm(direction_r);

                % Handle numerical instability
                if any(direction_r == Inf | direction_r == -Inf | ...
                       isnan(direction_r) | direction_r == -NaN)

                    r_3D = repmat(total_centroid_TMP(rogue_num_orig(end)+1,1:2), ...
                                  [numel(rogue_num_orig),1]); 
                    Interpolant_Result(blocks{1},1:2) = r_3D;

                else
                    % ------------------------------------------------------
                    % Project rogue points along fitted direction
                    % ------------------------------------------------------
                    midpoints = mean(rogue_volume);
                    dirvec    = direction_r';

                    dl = (rogue_num_orig - midpoints(end)) / dirvec(end);
                    dl = num2cell(dl,2);

                    points_slice = cellfun(@(x) ...
                        ([midpoints(1) dirvec(1); midpoints(2) dirvec(2)]*[1;x])', ...
                        dl,'uni',0);

                    points_slice = vertcat(points_slice{:});
                    r_3D = points_slice;

                    % Align extrapolated points with known centroid
                    Interpolant_Result(blocks{1},1:2) = ...
                        r_3D - repmat((r_3D(end,:) - ...
                        total_centroid_TMP(rogue_num_orig(end)+1,1:2)), ...
                        size(r_3D,1),1);
                end
            end
        end

        %% =================================================================
        % CASE 2: Rogue region at the END of the fiber
        % =================================================================
        if ismember(numel(non_rogue_area_find),find(non_rogue_area_find == 0))

            rogue_volume   = total_centroid_TMP(blocks{end},:);
            rogue_num_orig = rogue_volume(:,3);

            % Remove outliers if needed
            if numel(rogue_num_orig) > 1
                [rogue_volume,~] = regstats_func(rogue_volume);
                rogue_num_curr   = rogue_volume(:,3);
            else 
                rogue_num_curr   = rogue_num_orig;
            end

            % Fallback case
            if (numel(rogue_num_curr) == 1 || isempty(rogue_num_curr))

                r_3D = repmat(total_centroid_TMP(rogue_num_orig(1)-1,1:2), ...
                              [numel(rogue_num_orig),1]); 
                Interpolant_Result(blocks{end},1:2) = r_3D;

            else
                % Estimate direction via SVD
                r0        = mean(rogue_volume);
                xyz_cent  = bsxfun(@minus,rogue_volume,r0);
                [~,~,V]   = svd(xyz_cent,0);
                direction_r = V(:,1);
                direction_r = direction_r ./ vecnorm(direction_r);

                midpoints = mean(rogue_volume);
                dirvec    = direction_r';

                % Handle numerical instability
                if any(direction_r == Inf | direction_r == -Inf | ...
                       isnan(direction_r) | direction_r == -NaN)

                    r_3D = repmat(total_centroid_TMP(rogue_num_orig(1)-1,1:2), ...
                                  [numel(rogue_num_orig),1]); 
                    Interpolant_Result(blocks{end},1:2) = r_3D;

                else
                    % Project forward
                    dl = (rogue_num_orig - midpoints(end)) / dirvec(end);
                    dl = num2cell(dl,2);

                    points_slice = cellfun(@(x) ...
                        ([midpoints(1) dirvec(1); midpoints(2) dirvec(2)]*[1;x])', ...
                        dl,'uni',0);

                    points_slice = vertcat(points_slice{:});
                    r_3D = points_slice;

                    % Align extrapolated points
                    Interpolant_Result(blocks{end},1:2) = ...
                        r_3D - repmat((r_3D(1,:) - ...
                        total_centroid_TMP(rogue_num_orig(1)-1,1:2)), ...
                        size(r_3D,1),1);
                end
            end
        end
   end

   %% ======================================================================
   % Clamp coordinates to image bounds
   % Note: centroid format is (col,row) while size_length is (row,col)
   %% ======================================================================
   Interpolant_Result(Interpolant_Result(:,1) > size_length(2),1) = size_length(2);
   Interpolant_Result(Interpolant_Result(:,2) > size_length(1),2) = size_length(1);
   Interpolant_Result(Interpolant_Result(:,1) < 0,1) = 0.1;
   Interpolant_Result(Interpolant_Result(:,2) < 0,2) = 0.1;

   % ----------------------------------------------------------------------
   % Output corrected centroids
   % ----------------------------------------------------------------------
   varargout{1} = Interpolant_Result(:,1:2);

end