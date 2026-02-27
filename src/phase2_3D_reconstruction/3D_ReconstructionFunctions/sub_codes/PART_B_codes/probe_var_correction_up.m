function [varargout] = probe_var_correction_up(varargin)

% PROBE_VAR_CORRECTION_UP
% -------------------------------------------------------------------------
% Adjusts / reconstructs the centroid path when scanning upward.
% The function checks whether the fitted results extend beyond the
% original centroid range in the upward z-direction and performs
% interpolation if necessary.
%
% INPUTS (via varargin)
%   1) Centroid_MID_cell_LUMP_ORIG_UP : cell array of original centroids
%   2) c_f                             : current fiber index
%   3) probe_var                       : probe index used for z reference
%   4) main_res                        : fitted/estimated centroid results
%   5) z_TMP                           : cell array of z-slice locations
%
% OUTPUT
%   varargout{1} : corrected centroid trajectory (NEW_CENTROID)
% -------------------------------------------------------------------------

Centroid_MID_cell_LUMP_ORIG_UP     =    varargin{1};
c_f                                =    varargin{2};
probe_var                          =    varargin{3};
main_res                           =    varargin{4};
z_TMP                              =    varargin{5};

% Check if the fitted results extend above the probe's starting slice.
% If so, we may need to reconstruct the centroid path.
if size(main_res(main_res(:,3) < z_TMP{probe_var}(1),:),1) > size(Centroid_MID_cell_LUMP_ORIG_UP{c_f},1)
    
    % Case 1: Original centroid has multiple points → perform interpolation
    if size(Centroid_MID_cell_LUMP_ORIG_UP{c_f},1) > 1
        
        % Build an extended centroid path between:
        %   (a) the current fiber centroid path, and
        %   (b) the probe fiber starting point.
        % This enables interpolation instead of extrapolation.
        orig_cent_ext = [Centroid_MID_cell_LUMP_ORIG_UP{c_f} ; ...
                         Centroid_MID_cell_LUMP_ORIG_UP{probe_var}(1,:)];
        
        % Upward scan flag
        query_path_f   =  2;
        
        % Interpolate to fill missing slices
        NEW_CENTROID   =  MINOR_INTERPOLATION_FUNC(orig_cent_ext,query_path_f);
        
    % Case 2: Only one centroid available → use fitted results directly
    elseif size(Centroid_MID_cell_LUMP_ORIG_UP{c_f},1) == 1
        
        NEW_CENTROID = main_res(main_res(:,3) < z_TMP{probe_var}(1),:);
    end
    
else
    % No correction needed → retain the original centroid path
    NEW_CENTROID = Centroid_MID_cell_LUMP_ORIG_UP{c_f};
end

% Return corrected centroid path
varargout{1} = NEW_CENTROID;