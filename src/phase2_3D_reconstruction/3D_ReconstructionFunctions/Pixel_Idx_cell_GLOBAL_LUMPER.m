% =========================================================================
% PURPOSE:
% This function resolves ("lumps") nested Pixel_Idx_cell_GLOBAL structures
% into a single-level cell array indexed by numel_count.
%
% The input data may contain multiple tiers of cell nesting due to previous
% splitting/merging operations. This routine:
%   • Traverses all nesting levels
%   • Converts scalar cells to numeric arrays when needed
%   • Places each fiber's pixel index into its correct global position
%
% NOTE:
% The original control flow is intentionally preserved.
% =========================================================================

function [varargout] = Pixel_Idx_cell_GLOBAL_LUMPER(varargin)

% -------------------------------------------------------------------------
% INPUTS
% -------------------------------------------------------------------------
Pixel_Idx_cell_GLOBAL   =  varargin{1};  % Nested pixel index cell structure
numel_count             =  varargin{2};  % Target global index mapping


% -------------------------------------------------------------------------
% Initialize output container
% -------------------------------------------------------------------------
Pixel_Idx_cell_GLOBAL_LUMP = cell(1);

% -------------------------------------------------------------------------
% Main traversal over fibers
% -------------------------------------------------------------------------
for pp = 1:numel(numel_count)  
    
    % =====================================================================
    % CASE 1: numel_count{pp} itself is a cell (multi-tier structure)
    % =====================================================================
    if iscell(numel_count{pp})       
        
       for kk = 1:numel(numel_count{pp})
           
           % --------------------------------------------------------------
           % Subcase: multiple indices grouped together
           % --------------------------------------------------------------
           if numel(numel_count{pp}{kk}) > 1 
              
              for mm = 1:numel(numel_count{pp}{kk})
                  
                  % Convert scalar cell to numeric if needed
                  if numel(Pixel_Idx_cell_GLOBAL{pp}{kk}{mm}) == 1
                      Pixel_Idx_cell_GLOBAL_LUMP{numel_count{pp}{kk}(mm)} = ...
                          cell2mat(Pixel_Idx_cell_GLOBAL{pp}{kk}{mm}); 
                  else 
                      Pixel_Idx_cell_GLOBAL_LUMP{numel_count{pp}{kk}(mm)} = ...
                          Pixel_Idx_cell_GLOBAL{pp}{kk}{mm};    
                  end
              end 
              
           % --------------------------------------------------------------
           % Subcase: single index mapping
           % --------------------------------------------------------------
           else     
               if numel(Pixel_Idx_cell_GLOBAL{pp}{kk}) == 1
                   Pixel_Idx_cell_GLOBAL_LUMP{numel_count{pp}{kk}} = ...
                       cell2mat(Pixel_Idx_cell_GLOBAL{pp}{kk});
               else 
                   Pixel_Idx_cell_GLOBAL_LUMP{numel_count{pp}{kk}} = ...
                       Pixel_Idx_cell_GLOBAL{pp}{kk};    
               end
           end   
       end      
       
    % =====================================================================
    % CASE 2: numel_count{pp} is numeric (single-tier structure)
    % =====================================================================
    else       
        
        % --------------------------------------------------------------
        % Subcase: multiple target indices
        % --------------------------------------------------------------
        if numel((numel_count{pp})) > 1
            
            for kk = 1:numel(numel_count{pp})
                
                % Convert scalar cell to numeric if needed
                if numel(Pixel_Idx_cell_GLOBAL{pp}{kk}) == 1
                    Pixel_Idx_cell_GLOBAL_LUMP{numel_count{pp}(kk)} = ...
                        cell2mat(Pixel_Idx_cell_GLOBAL{pp}{kk}); 
                else 
                    Pixel_Idx_cell_GLOBAL_LUMP{numel_count{pp}(kk)} = ...
                        Pixel_Idx_cell_GLOBAL{pp}{kk};
                end
            end
            
        % --------------------------------------------------------------
        % Subcase: direct one-to-one mapping
        % --------------------------------------------------------------
        else 
            Pixel_Idx_cell_GLOBAL_LUMP{numel_count{pp}} = ...
                Pixel_Idx_cell_GLOBAL{pp};    
        end   
    end
end

% -------------------------------------------------------------------------
% OUTPUT
% -------------------------------------------------------------------------
varargout{1} = Pixel_Idx_cell_GLOBAL_LUMP;

end