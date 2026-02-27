
%%
clearvars; close all; clc;

% Load the .mat file 
img_vol_chunks_directory  = '/Users/ronaldagyei/Desktop/PhD_Codes/results/SPLIT_FILES';
mat_file_directory        = '/Users/ronaldagyei/Desktop/PhD_Codes/mat_files';


SPLIT_PATH_var = sprintf('%s%s',img_vol_chunks_directory,'/XYZ_files');

Input = load(sprintf('%s%s%s',mat_file_directory,'/','New_Grey_Data_S2.mat'));
Input = cell2mat(struct2cell(Input));

% Load the slice_stats_structure 
load(sprintf('%s%s%s',img_vol_chunks_directory,'/','FULL_VOLUME_METADATA_LIBRARY_0002_XYZ.mat'));

% tt = 100;
% imshow(Input(:,:,tt))
% SLICE_REGIONS_STATS_TEMP = SLICE_REGIONS_STATS_STRUCTURE{tt};

%%
% plot the ellipse 
tt              =   12;
thresh_max      =   6;
thresh_min      =   3;
Stats_FINAL     =   table2struct(SLICE_REGIONS_STATS_STRUCTURE{tt});
Image           =   Input(:,:,tt);
[xe,ye]         =   ellipse_plot_func_circ_post_process(Stats_FINAL,Image,tt,thresh_max,thresh_min);

%%





% Stats_check_Pixel  =    {Stats_FINAL.Pixel_IDX_List};
% Stats_check_Minor  =    {Stats_FINAL.MinorAxisLength};
% 
% %% Grey scale removal and MinorAxisLength removal 
% % pix_temp           =    cell2mat(cellfun(@(x)(numel(find(Image(x)>40000))/numel(x))>0.2,Stats_check_Pixel,'uni',0));
% 
% pix_temp           =    cell2mat(cellfun(@(x)(numel(find(Image(x)>40000))/numel(x))>0.2,Stats_check_Pixel,'uni',0));
% pix_temp           =    ~logical(pix_temp);
% 
% % Size removal 
% minor_temp         =    cell2mat(cellfun(@(x)x>15,Stats_check_Minor,'uni',0));
% minor_temp         =    logical(minor_temp);
% 
% Stats_FINAL(pix_temp|minor_temp)   =   [];
% [xe,ye]                            =   ellipse_plot_func_circ(Stats_FINAL,Image,73);
% 
% %% Employ if MinorAxisLength removal 
% 
% if isempty(Stats_FINAL)
% 
%     Stats_FINAL     =   SLICE_REGIONS_STATS_STRUCTURE{tt};    
%     Stats_check_Minor  =    cell2mat({Stats_FINAL.MinorAxisLength});       
%     Stats_FINAL(Stats_check_Minor > 15) = [];
% 
% end
% 
% %%
% tt              =   1;
% Image           =   Input_73_30(:,:,tt);
% 
% 
% 
% 
