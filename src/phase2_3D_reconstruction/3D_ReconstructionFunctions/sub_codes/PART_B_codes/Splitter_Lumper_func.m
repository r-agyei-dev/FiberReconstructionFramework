function [varargout] = Splitter_Lumper_func(Linear_Index_center_GLOBAL)

% SPLITTER_LUMPER_FUNC
% This function takes a possibly nested cell array of linear indices representing fibers
% and "flattens" it into a single-level cell array, where each entry corresponds
% to a single fiber segment. It handles:
%   - Nested cell arrays within TMP
%   - Empty sets
%   - Single-level arrays
% The result is a cleaned, lumped version of the fiber indices.

% Copy of the original input for processing
Linear_Index_center_GLOBAL_LUMP = Linear_Index_center_GLOBAL;

% Initialize output cell
Linear_Index_center_GLOBAL_LUMP_cell = cell(1,1);

% Initialize counter for output cell indexing
counter = 1;

% Loop through all fibers in the original data
for ii = 1:size(Linear_Index_center_GLOBAL,2)
    
    TMP = Linear_Index_center_GLOBAL{ii};  % Current fiber or group of fibers

    if iscell(TMP)
        % Check if TMP contains nested cells
        if any(cell2mat(cellfun(@(x)iscell(x),TMP,'uni',0)))
            TMP = [TMP{:}];  % Flatten nested cells
        end
        
        if iscell(TMP)
            % If still a cell array, iterate and add each sub-cell to output
            for k = 1:size(TMP,2)
                Linear_Index_center_GLOBAL_LUMP_cell{counter} = TMP{k};
                counter = counter + 1;                         
            end
        else
            % TMP is now a numeric array, add directly to output
            Linear_Index_center_GLOBAL_LUMP_cell{counter} = TMP;
        end
        
    else
        % TMP was not a cell, add directly to output
        Linear_Index_center_GLOBAL_LUMP_cell{counter} = TMP;
    end
    
    % Increment counter for next fiber
    counter = counter + 1;
    
end

% Assign the flattened lumped cell array to output
varargout{1} = Linear_Index_center_GLOBAL_LUMP_cell;


% -------------------------
% Notes / Reference:
% The commented section below shows an older version of the function, which
% attempted to flatten nested cells but did not handle all cases (empty cells,
% multi-level nesting). The current version safely flattens fibers for further
% processing.
%
% Example usage:
%   Lumped_Fibers = Splitter_Lumper_func(Linear_Index_center_GLOBAL);
%   Each entry Lumped_Fibers{i} now contains a single fiber segment's linear indices.