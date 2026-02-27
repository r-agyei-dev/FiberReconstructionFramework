function [varargout] = Rogue_region_FUNC(varargin)
% Rogue_region_FUNC
% -------------------------------------------------------------------------
% Identifies "rogue" regions based on the area distribution of elliptical regions.
% Returns a logical array indicating non-rogue regions (true) and rogue regions (false)
% based on comparison with adaptive thresholds derived from the mean size of 
% regions either above or below a set threshold.
%
% INPUTS (via varargin)
%   1  total_area_temp : array of areas of all candidate regions
%
% OUTPUTS
%   non_rogue_area_find : logical array (1 = non-rogue, 0 = rogue)
% -------------------------------------------------------------------------

total_area_temp = varargin{1};

% Define maximum area threshold for ellipses
max_ellip_thresh = 10;

% -------------------------------------------------------------------------
% Case 1: Some regions have area smaller than threshold
% -------------------------------------------------------------------------
if any(total_area_temp < max_ellip_thresh)
    
    tmp_var = 1.5; % multiplier for determining non-rogue areas
    
    % Extract areas smaller than threshold
    total_area_temp_extract = total_area_temp(total_area_temp < max_ellip_thresh);
    
    % Compute mean of values above the mean to focus on larger-than-average small regions
    mean_value_tmp = mean(total_area_temp_extract(total_area_temp_extract > mean(total_area_temp_extract)));
    
    % Determine non-rogue areas based on scaled mean
    non_rogue_area_find = round(total_area_temp) < round(tmp_var * mean_value_tmp,1);

% -------------------------------------------------------------------------
% Case 2: All regions exceed the threshold
% -------------------------------------------------------------------------
else
    tmp_var = 2; % multiplier for determining non-rogue areas
    
    % Extract areas greater than threshold
    total_area_temp_extract = total_area_temp(total_area_temp > max_ellip_thresh);
    
    % Compute mean of values below the mean to focus on smaller-than-average large regions
    mean_value_tmp = mean(total_area_temp_extract(total_area_temp_extract < mean(total_area_temp_extract)));
    
    % Determine non-rogue areas based on scaled mean
    non_rogue_area_find = round(total_area_temp) < round(tmp_var * mean_value_tmp,1);  
end

varargout{1} = non_rogue_area_find;

end