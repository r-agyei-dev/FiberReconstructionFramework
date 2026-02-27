function [SLICE_REGIONS_STATS_STRUCTURE] = METADATA_CROP_FUNC(SLICE_REGIONS_STATS_STRUCTURE,load_variable,crop_region)

% =========================================================================
% PURPOSE:
% Crop region metadata by removing objects whose centroids fall outside
% the specified crop bounds along a chosen axis.
%
% INPUTS:
%   SLICE_REGIONS_STATS_STRUCTURE : cell array of regionprops metadata
%   load_variable                 : axis selector (2 → xz, 3 → zy)
%   crop_region                   : [min max] bounds for cropping
%
% OUTPUT:
%   SLICE_REGIONS_STATS_STRUCTURE : updated structure with cropped regions
% =========================================================================

disp('Extracting the centroid from the cleaned tructure slice')

% -------------------------------------------------------------------------
% Extract centroids from each slice into a working cell array
% -------------------------------------------------------------------------
Centroid_cell = cell(1,length(SLICE_REGIONS_STATS_STRUCTURE));

for ii  =  1  :  length(SLICE_REGIONS_STATS_STRUCTURE)-1
    % Direct centroid extraction from structure
%     Centroid_cell{ii} = cat(1,(SLICE_REGIONS_STATS_STRUCTURE{ii}(:).Centroid));
    Centroid_cell{ii} = SLICE_REGIONS_STATS_STRUCTURE{ii}.Centroid;
end


% =========================================================================
% Perform cropping depending on the dataset orientation
% =========================================================================

if load_variable == 2
    
    % ---------------------------------------------------------------------
    % Case: XZ orientation
    % Centroid format: (cols, rows) → (x,z)
    % ---------------------------------------------------------------------
    for ii = 1:size(Centroid_cell,2)-1   
        
        % Identify centroids outside crop bounds (z-direction)
        del_idx = (Centroid_cell{ii}(:,2) < crop_region(1) | ...
                   Centroid_cell{ii}(:,2) > crop_region(2));

        % Remove corresponding regions from metadata
        SLICE_REGIONS_STATS_STRUCTURE{ii}(del_idx,:) = [];

        % Update centroid list
        Centroid_cell{ii}(del_idx,:) = [];
    end

elseif load_variable == 3
    
    % ---------------------------------------------------------------------
    % Case: ZY orientation
    % Centroid format: (cols, rows) → (z,y)
    % ---------------------------------------------------------------------
    for ii = 1:size(Centroid_cell,2)-1   
        
        % Identify centroids outside crop bounds (y-direction)
        del_idx = (Centroid_cell{ii}(:,1) < crop_region(1) | ...
                   Centroid_cell{ii}(:,1) > crop_region(2));

        % Remove corresponding regions from metadata
        SLICE_REGIONS_STATS_STRUCTURE{ii}(del_idx,:) = [];

        % Update centroid list
        Centroid_cell{ii}(del_idx,:) = [];
    end      
end



%%
% -------------------------------------------------------------------------
% Legacy parallel version (kept for reference)
% -------------------------------------------------------------------------
% function [SLICE_REGIONS_STATS_STRUCTURE,Centroid_cell] = METADATA_CROP_FUNC(Centroid_cell,SLICE_REGIONS_STATS_STRUCTURE,load_variable,crop_region)
%
% disp('Extracting the centroid from the cleaned tructure slice')
% Centroid_cell = cell(1,length(SLICE_REGIONS_STATS_STRUCTURE));
% parfor ii  =  1  :  length(SLICE_REGIONS_STATS_STRUCTURE)-1
%     Centroid_cell{ii} = cat(1,(SLICE_REGIONS_STATS_STRUCTURE{ii}(:).Centroid));
% end
%
% if load_variable == 2
%     for ii = 1:size(Centroid_cell,2)-1   
%         % Note Centroid_cell is in cols/rows (xz)    zxy is rows/cols/depth 
%         del_idx = (Centroid_cell{ii}(:,2)<crop_region(1) | Centroid_cell{ii}(:,2)>crop_region(2));
%         SLICE_REGIONS_STATS_STRUCTURE{ii}(del_idx)=[];
%         Centroid_cell{ii}(del_idx,:) = [];
%     end
% elseif load_variable == 3
%     for ii = 1:size(Centroid_cell,2)-1   
%         % Note Centroid_cell is in cols/rows (zy)    yzx is rows/cols/depth
%         del_idx = (Centroid_cell{ii}(:,1)<crop_region(1) | Centroid_cell{ii}(:,1)>crop_region(2));
%         SLICE_REGIONS_STATS_STRUCTURE{ii}(del_idx)=[];
%         Centroid_cell{ii}(del_idx,:) = [];
%     end      
% end