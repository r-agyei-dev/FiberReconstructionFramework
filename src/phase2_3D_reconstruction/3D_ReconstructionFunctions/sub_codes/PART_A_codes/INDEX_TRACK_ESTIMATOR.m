% THIS FUNCTION TRACKS THE HISTORY OF CENTROIDS AND STORES THEM AS A CELL ARRAY
% It monitors the stacking of ellipses across slices and updates the remaining ellipses to be tracked.
% 
% Sub-functions used:
% (a) sieve_track_func       - Tracks initial indices of ellipses in the first slice
% (b) sieve_track_func_more  - Tracks ellipses across subsequent slices

% Overview of steps:
% 1. Use sieve_track_func to track ellipses in the first slice. If the first index in
%    sieve_parameter_unique_CELL is zero, tracking stops at the current loop index. 
%    Otherwise, the function tracks until the first zero is encountered.
% 2. Update remnant_ellipses as ellipses are captured in each iteration.

function [track_MAIN_CELL, remnant_ellipses] = INDEX_TRACK_ESTIMATOR(varargin)

% Input arguments:
% numberOfImageFiles           = Total number of images (slices) in the 3D volume
% remnant_ellipses             = Array storing ellipses that remain to be stacked; updated dynamically
% sieve_parameter_unique_CELL  = Cell array of indices of ellipses already stacked
% iter_limit                   = Maximum number of slices to process
% numtrack                     = One less than the maximum number of slices, used for tracking in sieve_track_func_more

numberOfImageFiles            = varargin{1};
remnant_ellipses              = varargin{2};
sieve_parameter_unique_CELL   = varargin{3};
iter_limit                    = varargin{4};
numtrack                      = varargin{5};

% Step 1: Track the initial ellipses in the first slice
track_cell_first              = sieve_track_func(numberOfImageFiles, sieve_parameter_unique_CELL); 

% Initialize the main cell array to store tracks for all iterations
track_MAIN_CELL               = cell(1, iter_limit);  

% Step 2: Track ellipses across subsequent slices
for iter_run = 2 : iter_limit
    % Track ellipses between successive slices using sieve_track_func_more
    % Inputs: numtrack, current slice, remaining ellipses, and indices of stacked ellipses
    track_cell_further        = sieve_track_func_more(numtrack, iter_run, remnant_ellipses, sieve_parameter_unique_CELL);  
    
    % Store tracked ellipses for the current slice in the main cell
    track_MAIN_CELL{iter_run} = track_cell_further;  
end

% Step 3: Add the first slice tracking to the main cell array
track_MAIN_CELL{1} = track_cell_first;  