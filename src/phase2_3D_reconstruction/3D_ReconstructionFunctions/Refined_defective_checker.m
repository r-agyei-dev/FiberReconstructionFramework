% =========================================================================
% PURPOSE:
% Identify defective or incomplete fiber entries based on missing pixel data.
%
% The function scans each element of Pixel_Idx_cell_GLOBAL and flags cases
% where expected pixel index data is empty while the corresponding center
% index entry is a cell array.
%
% OUTPUT:
%   d_f(ii) = 1  → defective / incomplete entry detected
%   d_f(ii) = 0  → entry appears valid
% =========================================================================

function [varargout] = Refined_defective_checker(varargin)

% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
Linear_Index_center_GLOBAL  = varargin{1};   % Cell array of fiber centers
Pixel_Idx_cell_GLOBAL       = varargin{2};   % Nested pixel index cell array

% -------------------------------------------------------------------------
% Initialize defect flag vector (same length as fiber list)
% -------------------------------------------------------------------------
d_f = zeros(1,length(Pixel_Idx_cell_GLOBAL));

% -------------------------------------------------------------------------
% Loop through each fiber entry
% -------------------------------------------------------------------------
for ii = 1:length(Pixel_Idx_cell_GLOBAL)

   % Progress display (string is generated but not printed)
   sprintf('%s%d','checking ',ii)

   % Only evaluate entries where the center index is itself a cell
   if iscell(Linear_Index_center_GLOBAL{ii})

      % Check whether ANY nested pixel index cell is empty
      % cellfun -> tests each subcell for emptiness
      % cell2mat -> converts logical cell output to numeric array
      % any(...) -> true if at least one empty element exists
      if any(cell2mat(cellfun(@(x) isempty(x), ...
                               [Pixel_Idx_cell_GLOBAL{ii}{:}], ...
                               'uni',0)))

          % Flag this fiber as defective
          d_f(ii) = 1;
      end

   end

end

% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
varargout{1} = d_f;   % Defect indicator vector