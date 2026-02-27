function [varargout] = probe_var_correction_down(varargin)

% PROBE_VAR_CORRECTION_DOWN
% -------------------------------------------------------------------------
% Adjusts / reconstructs the centroid path when scanning downward.
% The function decides whether interpolation is required based on whether
% the downstream results (main_res) contain more slices than the original
% centroid track for the current fiber.
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% PRESERVE ORIGINAL CENTROID TRACK WHEN POSSIBLE %%%%%%%
% Compare the number of downstream fitted points with the number of
% original centroid points for the current fiber. If the fitted path
% extends further in z, interpolation may be required.

if size(main_res(main_res(:,3) > z_TMP{probe_var}(end),:),1) > size(Centroid_MID_cell_LUMP_ORIG_UP{c_f},1)
    
    % Case 1: Original centroid has multiple points → perform interpolation
    if size(Centroid_MID_cell_LUMP_ORIG_UP{c_f},1) > 1
        
        % Build extended centroid path using the probe reference endpoint
        orig_cent_ext = [Centroid_MID_cell_LUMP_ORIG_UP{probe_var}(end,:) ; ...
                         Centroid_MID_cell_LUMP_ORIG_UP{c_f}];
        
        % Downward scan flag
        query_path_b = 1;
        
        % Interpolate to fill missing slices
        NEW_CENTROID = MINOR_INTERPOLATION_FUNC(orig_cent_ext,query_path_b);
        
    % Case 2: Only a single original centroid → use fitted results directly
    elseif size(Centroid_MID_cell_LUMP_ORIG_UP{c_f},1) == 1
        NEW_CENTROID = main_res(main_res(:,3) > z_TMP{probe_var}(end),:);
    end
    
else
    % No extension needed → keep the original centroid path
    NEW_CENTROID = Centroid_MID_cell_LUMP_ORIG_UP{c_f};
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Return corrected centroid path
varargout{   1} = NEW_CENTROID;