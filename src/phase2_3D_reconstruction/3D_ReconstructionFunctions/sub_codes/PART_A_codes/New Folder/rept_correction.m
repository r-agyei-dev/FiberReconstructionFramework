
%% BACKGROUND TO THE CODE: 
% SCAN THROUGH THE REGIONS TO FIND THE VARIABLE AND REMOVE THE REPEATING TERMS 
% Here we scan through the terms of the structure file to find out the
% indices of the struture with overlapping centroids. This bug was probably caused from the centroid cell locate.

function[rept_cell_global,SLICE_UPDATE,CENTROID_UPDATE]   =   rept_correction(Centroid_cell,SLICE_UPDATE,CENTROID_UPDATE)

% Initialize the rept_cell_global variable 
rept_cell_global = cell(1,size(Centroid_cell,2));

for ii     =   1  :  size(Centroid_cell,2)
    clearvars rept_cell

pdist_ans         =  pdist2(Centroid_cell{ii},Centroid_cell{ii});

% Create a cell along the row of the pdist2 cell variable 
pdist_ans_cell    =   num2cell(pdist_ans,2);   

% Find the zeros of the withi each of these cells 
rept_cell         =   cellfun(@(x)find(x==0),pdist_ans_cell,'uni',0);

% for cells with more than one zeros, extract the 2nd to last terms of the cells 
rept_cell         =   cellfun(@(x)x(2:end),rept_cell,'uni',0);

% find the non-empty cells, here these are the ones with indices which are repetitive
% find the unique values of these indices and delete them from the structure file 

non_empty         =   cell2mat(cellfun(@(x)~isempty(x),rept_cell,'uni',0))>0;
if isempty(non_empty)
rept_cell         =   [];
else
rept_cell         =   unique(horzcat(rept_cell{non_empty}));

SLICE_UPDATE{ii}(rept_cell)        =  [];
CENTROID_UPDATE{ii}(rept_cell,:)   =  [];
end
end



% This is a verification check part of the script.
% concantenate everythin into one 
non_empty = cell2mat(cellfun(@(x)~isempty(x),rept_cell_global,'uni',0))>0;
if isempty(non_empty)
rept_cell_global = [];
else
rept_cell_global = unique(horzcat(rept_cell_global{non_empty}));

end

