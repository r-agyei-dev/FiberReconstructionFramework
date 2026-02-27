%% ALIGNMENT_SUB_FUNCTION_CODE
% This function organizes and separates fiber segments in regions where fibers split
% (rogue regions) and stores their centroid coordinates for further processing.

function [varargout] = ALIGNMENT_SUB_FUNCTION_CODE(varargin)

% INPUTS:
% ALIGNMENT_SPLIT_FIBERS -> matrix of start and end indices of fiber split regions
% total_centroid_temp    -> centroid coordinates for all pixels
% Final_Result           -> main fiber coordinates after initial processing
% r_3D                   -> 3D extrapolated points for fiber splits
% blocks                 -> auxiliary fiber segment blocks (not directly used here)

% OUTPUTS:
% varargout{1} -> cell array of pixel coordinates for each distinct fiber segment
% varargout{2} -> concatenated pixels for all segments (PIX_TOT)
% varargout{3} -> dimensions (number of pixels) for each segment (pix_dim)

ALIGNMENT_SPLIT_FIBERS = varargin{1};
total_centroid_temp    = varargin{2};   
Final_Result           = varargin{3};  
r_3D                   = varargin{4};
blocks                 = varargin{5};

if ~isempty(ALIGNMENT_SPLIT_FIBERS)
    
    if numel(blocks) > 1  % Check if multiple fiber blocks exist

        % Adjust indices to include rogue regions at start and end
        ALIGNMENT_SPLIT_FIBERS = ALIGNMENT_SPLIT_FIBERS + repmat([1 -1], size(ALIGNMENT_SPLIT_FIBERS,1),1);

        % Convert split fiber indices into cell arrays for easier processing
        ALIGNMENT_SPLIT_FIBERS_cell = cell(1, size(ALIGNMENT_SPLIT_FIBERS,1));
        for ii = 1:size(ALIGNMENT_SPLIT_FIBERS,1)
            ALIGNMENT_SPLIT_FIBERS_cell{ii} = ALIGNMENT_SPLIT_FIBERS(ii,1) : ALIGNMENT_SPLIT_FIBERS(ii,2);
        end

        % Zero out shared regions in the main fiber matrix
        Final_Result([ALIGNMENT_SPLIT_FIBERS_cell{:}], :) = 0;

        % Create a temporary indexing matrix for splitting pixels
        clearvars ex_temp
        ex_temp(1,:) = ones(1, length(total_centroid_temp));
        ex_temp(2,:) = 1:length(total_centroid_temp);
        ex_temp(1,[ALIGNMENT_SPLIT_FIBERS_cell{:}]) = 0;

        % Split continuous non-rogue segments into separate cells
        num_z = sum(ex_temp(1,:) == 0);
        ex_temp_cell = cell(1,num_z);
        for tt = 1:num_z
            tmp_zero_idx = find(ex_temp(1,:) == 0, 1, 'first');
            ex_temp_cell{tt} = ex_temp(end, 1:tmp_zero_idx-1);
            ex_temp(:, 1:tmp_zero_idx) = [];
        end
        ex_temp_cell{tt+1} = ex_temp(2,:);
        ex_temp_cell = ex_temp_cell(~cellfun('isempty', ex_temp_cell));

        % Compute the length of each rogue region to be extrapolated
        length_to_extrapolate = ALIGNMENT_SPLIT_FIBERS(:,2) - ALIGNMENT_SPLIT_FIBERS(:,1) + 1;

        % Preallocate storage for pixels, concatenated pixels, and segment dimensions
        pixels  = cell(1, length(length_to_extrapolate)+1);
        PIX_TOT = cell(1, length(length_to_extrapolate)+1);
        pix_dim = cell(1, length(length_to_extrapolate)+1);

        % Loop through all fiber segments (including rogue splits)
        for ii = 1:length(length_to_extrapolate)+1

            % Store centroid coordinates for each segment along with vestiges from splits
            if ii == 1
                pixels{ii}{1} = Final_Result(ex_temp_cell{ii}, :);
                pixels{ii}{2} = r_3D{ii}{1};
                PIX_TOT{ii}   = vertcat(pixels{ii}{:});
                pix_dim{ii}   = cellfun(@(x) size(x,1), pixels{ii});

            elseif ii > 1 && ii < length(length_to_extrapolate)+1
                pixels{ii}{1} = r_3D{ii-1}{2};
                pixels{ii}{2} = Final_Result(ex_temp_cell{ii}, :);
                pixels{ii}{3} = r_3D{ii}{1};
                PIX_TOT{ii}   = vertcat(pixels{ii}{:});
                pix_dim{ii}   = cellfun(@(x) size(x,1), pixels{ii});

            elseif ii == length(length_to_extrapolate)+1
                pixels{ii}{1} = r_3D{ii-1}{2};
                pixels{ii}{2} = Final_Result(ex_temp_cell{ii}, :);
                PIX_TOT{ii}   = vertcat(pixels{ii}{:});
                pix_dim{ii}   = cellfun(@(x) size(x,1), pixels{ii});
            end
        end

    else
        % Case when only one fiber block exists (no splits)
        pixels  = total_centroid_temp;
        PIX_TOT = pixels;
        pix_dim = length(Final_Result);
    end

else
    % If no splits exist, return the original fiber result
    pixels  = Final_Result;
    PIX_TOT = Final_Result;
    pix_dim = length(Final_Result);
end

% Return outputs
varargout{1} = pixels;
varargout{2} = PIX_TOT;
varargout{3} = pix_dim;

end