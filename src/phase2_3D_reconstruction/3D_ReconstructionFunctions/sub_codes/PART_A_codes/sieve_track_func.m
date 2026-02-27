% THIS FUNCTION TRACKS ELLIPSE INDICES ACROSS SLICES STARTING FROM THE FIRST SLICE
%
% Purpose:
%   - For each ellipse in the first slice, track its successive indices through the 3D volume
%   - If the starting index in sieve_parameter_unique_CELL is zero, tracking stops immediately
%   - Otherwise, continue tracking along successive slices until a zero index is encountered
%   - Stores the tracked sequence for each ellipse in a cell array

function [track_cell] = sieve_track_func(numberOfImageFiles, sieve_parameter_unique_CELL)

% Initialize the cell array to store the tracked sequences for each ellipse in the first slice
track_cell = cell(1, length(sieve_parameter_unique_CELL{1}));

% Loop over each ellipse in the first slice
for tt = 1 : size(sieve_parameter_unique_CELL{1}, 1)

    % Start with the current ellipse index
    track_temp = tt;
    
    % Initialize the tracking sequence with the first index
    track(1) = track_temp;

    % Case 1: If the first index is zero, no further tracking is needed
    if sieve_parameter_unique_CELL{1}(track(1)) == 0
        track_cell{tt} = tt;

    % Case 2: If the first index is non-zero, track successive indices until zero is reached
    else
        for ii = 1 : numberOfImageFiles-1
            track_temp = sieve_parameter_unique_CELL{ii}(track_temp, end);  % Follow the next index
            track(ii+1) = track_temp;                                       % Append to the tracking sequence
            if track_temp == 0
                break                                                        % Stop tracking when zero encountered
            end
        end

        % Save the complete tracked sequence for this ellipse
        track_cell{tt} = track;

    end

    % Clean up temporary variables
    clearvars track track_temp
end