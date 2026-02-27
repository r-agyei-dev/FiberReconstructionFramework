function [varargout] = Rogue_region_FUNC_2(varargin)
% Rogue_region_FUNC_2
% -------------------------------------------------------------------------
% Identifies "rogue" regions based on area thresholds of elliptical regions.
% The function returns a logical array indicating non-rogue areas (true) and
% rogue areas (false). Thresholds are determined adaptively based on mean 
% and standard deviation of the area distribution.
%
% INPUTS (via varargin)
%   1  total_area_temp : array of areas of all candidate regions
%   2  strict_tmp_val  : strict ratio threshold for deciding rogue regions
%   3  strict_var      : strict multiplier for standard deviation
%   4  ease_var        : relaxed multiplier for standard deviation
%
% OUTPUTS
%   non_rogue_area_find : logical array (1 = non-rogue, 0 = rogue)
% -------------------------------------------------------------------------

total_area_temp    = varargin{1};
strict_tmp_val     = varargin{2};
strict_var         = varargin{3};
ease_var           = varargin{4};

% Define maximum area threshold for an ellipse
max_ellip_thresh = 10;

% -------------------------------------------------------------------------
% If there are ellipses exceeding the max threshold, apply conditional logic
% to determine which regions are considered non-rogue
% -------------------------------------------------------------------------
if max(total_area_temp) > max_ellip_thresh
    
    % (a) Small fraction of regions exceed threshold OR
    % (b) Large fraction exceed threshold but very few are below threshold
    if ((nnz(total_area_temp > max_ellip_thresh)/numel(total_area_temp)) <= strict_tmp_val) || ...
       ((nnz(total_area_temp > max_ellip_thresh)/numel(total_area_temp)) > strict_tmp_val && nnz(total_area_temp < max_ellip_thresh) <= 1)
   
        % Non-rogue regions are those below mean + strict multiplier * std
        non_rogue_area_find = total_area_temp < round(mean(total_area_temp) + strict_var*std(total_area_temp),1); 
        
    % (c) Large fraction exceed threshold and multiple regions below threshold
    elseif (nnz(total_area_temp > max_ellip_thresh)/numel(total_area_temp)) > strict_tmp_val && nnz(total_area_temp < max_ellip_thresh) > 1
        
        % Closely ranged rogue regions: max area < 2 * threshold
        if max(total_area_temp) < 2*max_ellip_thresh
            non_rogue_area_find = total_area_temp < max_ellip_thresh;
        
        % Sparsely ranged rogue regions: use relaxed threshold
        else
            non_rogue_area_find = total_area_temp < round(mean(total_area_temp) + ease_var*std(total_area_temp),1);
        end
        
    end

% If no ellipse exceeds the threshold, apply relaxed threshold
else
    non_rogue_area_find = total_area_temp < round(mean(total_area_temp) + ease_var*std(total_area_temp),1);
end

varargout{1} = non_rogue_area_find;

end

% -------------------------------------------------------------------------
% Notes / example usage:
% below_extract = total_area_temp(total_area_temp<10);
% mean(below_extract(below_extract>mean(below_extract)))
% This can be used to examine distribution of small ellipses for debugging.
% -------------------------------------------------------------------------