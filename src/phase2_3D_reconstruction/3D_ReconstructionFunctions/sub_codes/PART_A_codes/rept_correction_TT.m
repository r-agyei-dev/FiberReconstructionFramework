%% BACKGROUND:
% THIS FUNCTION SCANS THROUGH THE FIBER REGIONS TO IDENTIFY AND REMOVE
% DUPLICATE OR OVERLAPPING CENTROIDS
% 
% Purpose:
%   - Loop through each slice of the 3D volume
%   - Identify overlapping pixels among ellipses
%   - Remove duplicate entries from SLICE_UPDATE and CENTROID_UPDATE
%   - Ensure each centroid corresponds to a unique fiber region

function [SLICE_UPDATE, CENTROID_UPDATE] = rept_correction_TT(Centroid_cell, SLICE_UPDATE, CENTROID_UPDATE, size_length)

% Loop through each slice in the volume
for ii = 1 : size(Centroid_cell, 2)

    % Initialize a 2D map of the current slice to track pixels occupied by ellipses
    TT = zeros([size_length(1), size_length(2)]);

    % Mark the pixels of each ellipse in TT
    for tt = 1 : size(SLICE_UPDATE{ii}, 1)
        TT(SLICE_UPDATE{ii}(tt).Pixel_IDX_List) = tt;  % Assign ellipse index to occupied pixels
    end

    % Find indices of ellipses that have no unique pixels (duplicates)
    del_idx = ~ismember(1:size(SLICE_UPDATE{ii}), unique(nonzeros(TT)));

    % Remove duplicate ellipses from SLICE_UPDATE and corresponding centroids
    SLICE_UPDATE{ii}(del_idx)      = [];
    CENTROID_UPDATE{ii}(del_idx, :) = [];
    
end

end

%% NOTES ON UNUSED/COMMENTED SECTIONS:
% The commented code below provides alternative verification methods:
%   - Checks for repeated centroid coordinates
%   - Uses pdist2 to compute pairwise distances between centroids
%   - Identifies duplicates by zero distances and removes them
%   - Maintains a global record of all repetitive indices in rept_cell_global
%
% These sections are left for reference or debugging and are not active in
% the main function.
    
    
%%    
%     yes_find = zeros(1,351);
%     for ii = 1:351
%  Cent_temp =  Centroid_cell{ii};
%  
% clearvars val_idx 
% [val_idx(:,2),val_idx(:,1)] =  hist(Cent_temp(:,1),unique(Cent_temp(:,1))); 
%     
% comp_1 = ismember(Cent_temp(:,1),val_idx(val_idx(:,2)>1,1));
% 
% clearvars val_idx
% [val_idx(:,2),val_idx(:,1)] =  hist(Cent_temp(:,2),unique(Cent_temp(:,2))); 
%     
% comp_2 = ismember(Cent_temp(:,2),val_idx(val_idx(:,2)>1,1));
% 
% 
% if any(comp_1 & comp_2)
%    yes_find(ii) = 1; 
% end
%     end
% 
% 
% 
% Centroid_cell{ii}(235,:) = Centroid_cell{ii}(101,:);
% Centroid_cell{ii}(236,:) = Centroid_cell{ii}(118,:);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
%     
% pdist_ans         =  pdist2(Centroid_cell{ii},Centroid_cell{ii});
% 
% % Create a cell along the row of the pdist2 cell variable 
% pdist_ans_cell    =   num2cell(pdist_ans,2);   
% 
% % Find the zeros of the withi each of these cells 
% rept_cell         =   cellfun(@(x)find(x==0),pdist_ans_cell,'uni',0);
% 
% % for cells with more than one zeros, extract the 2nd to last terms of the cells 
% rept_cell         =   cellfun(@(x)x(2:end),rept_cell,'uni',0);
% 
% % find the non-empty cells, here these are the ones with indices which are repetitive
% % find the unique values of these indices and delete them from the structure file 
% 
% non_empty         =   cell2mat(cellfun(@(x)~isempty(x),rept_cell,'uni',0))>0;
% if isempty(non_empty)
% rept_cell         =   [];
% else
% rept_cell         =   unique(horzcat(rept_cell{non_empty}));
% 
% SLICE_UPDATE{ii}(rept_cell)        =  [];
% CENTROID_UPDATE{ii}(rept_cell,:)   =  [];
% end
% end
% 
% 
% 
% % This is a verification check part of the script.
% % concantenate everythin into one 
% non_empty = cell2mat(cellfun(@(x)~isempty(x),rept_cell_global,'uni',0))>0;
% if isempty(non_empty)
% rept_cell_global = [];
% else
% rept_cell_global = unique(horzcat(rept_cell_global{non_empty}));
% 
% end
% 
