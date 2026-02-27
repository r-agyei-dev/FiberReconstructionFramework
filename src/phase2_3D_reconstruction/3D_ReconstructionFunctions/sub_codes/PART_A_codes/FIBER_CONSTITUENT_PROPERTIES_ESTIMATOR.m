% THIS FUNCTION CREATES A CELL ARRAY OF METADATA FOR THE CONSTITUENT ELLIPSES OF FIBERS
% It consolidates all metadata of candidate fiber regions and reorders them 
% based on their spatial positions along consecutive slices.

function [properties] = FIBER_CONSTITUENT_PROPERTIES_ESTIMATOR (varargin)

% Input arguments:
% track_MAIN_CELL               = Cell array containing indices of ellipses stacked across consecutive slices for each fiber
% SLICE_REGIONS_STATS_STRUCTURE = Structure array containing metadata (properties) of all ellipses per slice

track_MAIN_CELL               = varargin{1};
SLICE_REGIONS_STATS_STRUCTURE = varargin{2};

% Initialize the output cell array to store the metadata of each fiber's constituent ellipses
properties = cell(1,length(track_MAIN_CELL));

% Loop through each fiber in track_MAIN_CELL
for hh = 1:length(track_MAIN_CELL)

    if hh == 1
        % First fiber: extract and group metadata directly from the slices
        for ii = 1:length(track_MAIN_CELL{hh})
            for kk = 1:length(track_MAIN_CELL{hh}{ii})  
                if track_MAIN_CELL{hh}{ii}(kk) == 0
                    break  % Stop if a zero marker is found
                else
                    % Assign the metadata of the kk-th ellipse in the ii-th group of the hh-th fiber
                    properties{hh}{ii}(kk,:) = SLICE_REGIONS_STATS_STRUCTURE{kk}(track_MAIN_CELL{hh}{ii}(kk),:);  
                end
            end
        end

    else
        % For subsequent fibers: adjust slice index to maintain spatial ordering
        for ii = 1:length(track_MAIN_CELL{hh})
            for kk = 1:length(track_MAIN_CELL{hh}{ii}) 
                if track_MAIN_CELL{hh}{ii}(kk) == 0
                    break  % Stop if a zero marker is found
                else
                    % Assign metadata while accounting for fiber's relative position in slice sequence
                    properties{hh}{ii}(kk,:) = SLICE_REGIONS_STATS_STRUCTURE{(hh+kk)-1}(track_MAIN_CELL{hh}{ii}(kk),:);  
                end
            end
        end
    end    
end