% THIS FUNCTION TRACKS THE HISTORY OF INDICES FOR ELLIPSES ACROSS SLICES
%
% Purpose:
%   - Track the sequence of ellipse indices starting from a given slice
%   - If the first index in sieve_parameter_unique_CELL is zero, tracking stops immediately
%   - Otherwise, continue tracking successive indices until a zero is encountered
%   - Similar to sieve_track_func but adapted for different indexing requirements

function [track_cell] = sieve_track_func_more(numberOfImageFiles, iter_run, remnant_ellipses, sieve_parameter_unique_CELL)

% Initialize the cell array to store tracked sequences for each ellipse in the current slice
track_cell = cell(1, length(remnant_ellipses{iter_run}));

% Loop over each ellipse in the current slice
for tt = 1 : length(remnant_ellipses{iter_run})
    
    % Start with the current ellipse index
    track_temp = remnant_ellipses{iter_run}(tt);

    % Initialize the tracking sequence with the first index
    track(1) = track_temp;

    % Case 1: If the first index is zero, no further tracking is needed
    if sieve_parameter_unique_CELL{iter_run}(track(1)) == 0
        track_cell{tt} = track_temp;

    % Case 2: If the first index is non-zero, track successive indices until zero is reached
    else
        for aa = iter_run : numberOfImageFiles-1
            track_temp = sieve_parameter_unique_CELL{aa}(track_temp, end);  % Follow the next index
            track(aa - iter_run + 2) = track_temp;                           % Append to the sequence
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