                                                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                            % This is a snippet of the code for the compariosn between 
                                                                            % two loading steps that uses Arithmetic operations only
                                                                            % BY: Ronald F Agyei.
                                                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear; close all; clc;
                                                                          
MUTUAL_nom_cell   =   {'MUTUAL_RESULT_OPTIMIZED_SEIVEDn'};
REFINED_nom_cell  =   {'REFINED_FIBERS_SIEVED'};
TEMP_CRUDE_cell   =   'TEMP_CRUDE_FULL_';
CENT_COORD_cell   =   {'Cent_Coord_UNSIEVED'};
load_steps_idx    =    73;
sieve_indicator   =    1;
pp                =    1;

for pp = 1
    Center_Idx_prefix         =    {'Center_Idx_SIEVED','Center_Idx_UNSIEVED'};
    initial_directory_name    =    sprintf('%s%d%s','/scratch/halstead/r/ragyei/CRITICAL_REGIONS_CELL/',load_steps_idx(pp),'_files');  % Only make changes to this path 

    %% initiate the while loop 
    clearvars MAIN_COUNT

    loop_track = [12 13 23];

    Idx_1   =   cell(1,3);
    Idx_2   =   cell(1,3);
    Idx_3   =   cell(1,3);
    
    for ii = 1:3

    track_seq = loop_track(ii);

    if track_seq == 12
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%         VOL_PROBE = REFINED_XYZ;     VOL_QUERY = REFINED_ZXY;      
        cd('/scratch/halstead/r/ragyei/Mutual_Comparison_GITHUB_LINUX')        
        [REFINED_XYZ,REFINED_ZXY] = H5_Loader_Function(initial_directory_name,load_steps_idx(pp),sieve_indicator,REFINED_nom_cell,track_seq);
        size_length = size(REFINED_XYZ);
        
        cd('/scratch/halstead/r/ragyei/Mutual_Comparison_GITHUB_LINUX')
%         tic
        disp('begin EPAD_FUNCTION_CODE_MASTER_VERSION_MUT')
        [MAIN_COUNT_probe12,MAIN_COUNT_query12,SINGLE_CELL_12] = EPAD_FUNCTION_CODE_MASTER_VERSION_MUT(REFINED_XYZ,REFINED_ZXY);
%         toc

%         disp('saving the variables')
%         save(fullfile(initial_directory_name,sprintf('MAIN_COUNT_probe12.mat')),'MAIN_COUNT_probe12','-v7.3');
%         save(fullfile(initial_directory_name,sprintf('MAIN_COUNT_query12.mat')),'MAIN_COUNT_query12','-v7.3');
%         save(fullfile(initial_directory_name,sprintf('SINGLE_CELL_12.mat')),'SINGLE_CELL_12','-v7.3');
        
        clearvars   REFINED_XYZ       REFINED_ZXY 
        [Idx_1{ii},Idx_2{ii}] = MUTUAL_COMPARISON_CODE_ANALYZERS(MAIN_COUNT_probe12,MAIN_COUNT_query12,SINGLE_CELL_12);
        

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        IDX_PROBE = LIDX_XYZ;      IDX_QUERY = LIDX_ZXY; slice_PROBE = slice_XYZ;      slice_QUERY = slice_ZXY;     
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    elseif track_seq == 13
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    
        [REFINED_XYZ,REFINED_YZX] = H5_Loader_Function(initial_directory_name,load_steps_idx(pp),sieve_indicator,REFINED_nom_cell,track_seq);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cd('/scratch/halstead/r/ragyei/Mutual_Comparison_GITHUB_LINUX')
%         tic
        disp('begin EPAD_FUNCTION_CODE_MASTER_VERSION_MUT')
        [MAIN_COUNT_probe13,MAIN_COUNT_query13,SINGLE_CELL_13] = EPAD_FUNCTION_CODE_MASTER_VERSION_MUT(REFINED_XYZ,REFINED_YZX);
%         toc

%         disp('saving the variables')
%         save(fullfile(initial_directory_name,sprintf('MAIN_COUNT_probe13.mat')),'MAIN_COUNT_probe13','-v7.3');
%         save(fullfile(initial_directory_name,sprintf('MAIN_COUNT_query13.mat')),'MAIN_COUNT_query13','-v7.3');
%         save(fullfile(initial_directory_name,sprintf('SINGLE_CELL_13.mat')),'SINGLE_CELL_13','-v7.3');
        
        clearvars REFINED_XYZ      REFINED_YZX
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        [Idx_1{ii},Idx_3{ii}] = MUTUAL_COMPARISON_CODE_ANALYZERS(MAIN_COUNT_probe13,MAIN_COUNT_query13,SINGLE_CELL_13);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    elseif track_seq == 23
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        cd('/scratch/halstead/r/ragyei/Mutual_Comparison_GITHUB_LINUX')
        [REFINED_ZXY,REFINED_YZX] = H5_Loader_Function(initial_directory_name,load_steps_idx(pp),sieve_indicator,REFINED_nom_cell,track_seq);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cd('/scratch/halstead/r/ragyei/Mutual_Comparison_GITHUB_LINUX')
%         tic
        disp('begin EPAD_FUNCTION_CODE_MASTER_VERSION_MUT')        
        [MAIN_COUNT_probe23,MAIN_COUNT_query23,SINGLE_CELL_23] = EPAD_FUNCTION_CODE_MASTER_VERSION_MUT(REFINED_ZXY,REFINED_YZX);
%         toc

%         disp('saving the variables')
%         save(fullfile(initial_directory_name,sprintf('MAIN_COUNT_probe23.mat')),'MAIN_COUNT_probe23','-v7.3');
%         save(fullfile(initial_directory_name,sprintf('MAIN_COUNT_query23.mat')),'MAIN_COUNT_query23','-v7.3');
%         save(fullfile(initial_directory_name,sprintf('SINGLE_CELL_23.mat')),'SINGLE_CELL_23','-v7.3');        
        clearvars REFINED_ZXY     REFINED_YZX       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        [Idx_2{ii},Idx_3{ii}] = MUTUAL_COMPARISON_CODE_ANALYZERS(MAIN_COUNT_probe23,MAIN_COUNT_query23,SINGLE_CELL_23);
                        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    end
    end

%   save(fullfile(initial_directory_name,sprintf('Idx_1.mat')),'Idx_1','-v7.3');
%   save(fullfile(initial_directory_name,sprintf('Idx_2.mat')),'Idx_2','-v7.3');
%   save(fullfile(initial_directory_name,sprintf('Idx_3.mat')),'Idx_3','-v7.3');
        
Idx_1_fin   =   [Idx_1{:}];
Idx_2_fin   =   [Idx_2{:}];
Idx_3_fin   =   [Idx_3{:}];

%% Redefine the Linear Index of the fibers for the 3D volume:
LIDX       =   load(sprintf('%s%d%s%d%s','/scratch/halstead/r/ragyei/CRITICAL_REGIONS_CELL/',load_steps_idx(pp),'_files/XYZ_Recon_files_CONDENSED_N_',load_steps_idx(pp),'/Center_Idx_UNSIEVED_xyz.mat'));
names      =   fieldnames(LIDX);
LIDX       =   LIDX.(names{1});
LIDX_var   =   LIDX(Idx_1_fin);

clearvars LIDX 
LIDX       =   load(sprintf('%s%d%s%d%s','/scratch/halstead/r/ragyei/CRITICAL_REGIONS_CELL/',load_steps_idx(pp),'_files/XYZ_Recon_files_CONDENSED_N_',load_steps_idx(pp),'/Cent_Coord_UNSIEVED_xyz.mat'));
names      =   fieldnames(LIDX);
LIDX       =   LIDX.(names{1});
Cent_var   =   LIDX(Idx_1_fin);

clearvars LIDX 
LIDX       =    load(sprintf('%s%d%s%d%s','/scratch/halstead/r/ragyei/CRITICAL_REGIONS_CELL/',load_steps_idx(pp),'_files/ZXY_Recon_files_CONDENSED_N_',load_steps_idx(pp),'/Center_Idx_UNSIEVED_zxy.mat'));
names      =    fieldnames(LIDX);
LIDX       =    LIDX.(names{1});
LIDX_var   =   [LIDX_var     LIDX(Idx_2_fin)];

clearvars LIDX 
LIDX       =    load(sprintf('%s%d%s%d%s','/scratch/halstead/r/ragyei/CRITICAL_REGIONS_CELL/',load_steps_idx(pp),'_files/ZXY_Recon_files_CONDENSED_N_',load_steps_idx(pp),'/Cent_Coord_UNSIEVED_zxy.mat'));
names      =    fieldnames(LIDX);
LIDX       =    LIDX.(names{1});
Cent_var   =   [Cent_var     LIDX(Idx_2_fin)];


clearvars LIDX
LIDX       =   load(sprintf('%s%d%s%d%s','/scratch/halstead/r/ragyei/CRITICAL_REGIONS_CELL/',load_steps_idx(pp),'_files/YZX_Recon_files_CONDENSED_N_',load_steps_idx(pp),'/Center_Idx_UNSIEVED_yzx.mat'));
names      =   fieldnames(LIDX);
LIDX       =   LIDX.(names{1});
LIDX_var   =   [LIDX_var     LIDX(Idx_3_fin)];

clearvars LIDX
LIDX       =   load(sprintf('%s%d%s%d%s','/scratch/halstead/r/ragyei/CRITICAL_REGIONS_CELL/',load_steps_idx(pp),'_files/YZX_Recon_files_CONDENSED_N_',load_steps_idx(pp),'/Cent_Coord_UNSIEVED_yzx.mat'));
names      =   fieldnames(LIDX);
LIDX       =   LIDX.(names{1});
Cent_var   =   [Cent_var     LIDX(Idx_3_fin)];



TEMP_MAT = zeros(size_length);

for hh = 1:length(LIDX_var) %142949 %142950%266655%length(LIDX_var)
    if all(TEMP_MAT(LIDX_var{hh}) == 0)
       if numel(LIDX_var{hh})>=100 
	   TEMP_MAT(LIDX_var{hh}) = hh;
       end
    else        
%   break
        unique_idx  =  unique(nonzeros(TEMP_MAT(LIDX_var{hh})));
        tmp_idx     =  vertcat(unique_idx,hh);
        LIDX_tmp    =  LIDX_var(tmp_idx);
%         [~,~,z]     =  cellfun(@(x)ind2sub(size_length,x),LIDX_tmp,'uni',0);
        z           =  cell2mat(cellfun(@(x)numel(x),LIDX_tmp,'uni',0));    
        [~,max_z] = max(z);
    
    if tmp_idx(max_z) == hh    
    TEMP_MAT(vertcat(LIDX_tmp{:})) = 0;
%     LIDX_tmp(~ismember(1:numel(z),max_z)) = [];
    TEMP_MAT(LIDX_tmp{max_z}) = hh;  
    end
        
%     % Update the variable LIDX_var
%     [LIDX_var{z(~ismember(1:numel(z),max_z))}] = deal({});
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% STEPS for ensuring accurate 3D fiber reconstructions 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% STEP 1: BINARIZE IMAGE VOLUME AND PERFORM SLICE-WISE REMOVAL OF ALL THE SMALL FEATURES ALONG ALL THE DIMENSIONS 
% BINARIZE 
TEMP_BIN = TEMP_MAT ~= 0;

% SMALL FEATURES REMOVAL ALONG THE XYZ DXN
for ii = 1:size(TEMP_BIN,3)    
TEMP_BIN(:,:,ii)  = bwareaopen(TEMP_BIN(:,:,ii),5);      
end

% SMALL FEATURES REMOVAL ALONG THE ZXY DXN
TEMP_BIN                =   permute(TEMP_BIN,[3 1 2]);              % (z x y)
for ii = 1:size(TEMP_BIN,3)    
TEMP_BIN(:,:,ii)  = bwareaopen(TEMP_BIN(:,:,ii),5);      
end
TEMP_BIN              =   permute(TEMP_BIN,[2 3 1]);                % (z x y) ==>> xyz
% 
% SMALL FEATURES REMOVAL ALONG THE YZX DXN
TEMP_BIN                =   permute(TEMP_BIN,[2 3 1]);              % (y z x)
for ii = 1:size(TEMP_BIN,3)    
TEMP_BIN(:,:,ii)  = bwareaopen(TEMP_BIN(:,:,ii),5);      
end
TEMP_BIN               =   permute(TEMP_BIN,[3 1 2]);               % (y z x) ==>> xyz
%

% REMOVE THE 3D USING THE SLICE REMOVAL MASK
TEMP_MAT = TEMP_MAT.*TEMP_BIN;

%% STEP 2: REMOVE ANY FIBER THAT MAY BE WITHIN THE BAND OF FALSE RECONSTRUCTION AROUND THE SPECIMEN 

imageSizeX = size(TEMP_MAT,2);  imageSizeY = size(TEMP_MAT,1);
[columnsInImage,rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% Next create the circle in the image.
centerX = 0.5*imageSizeX;  centerY = 0.5*imageSizeY;

% if pp == 1
    edge_radius_2 = 0.97*centerX;
    edge_removal = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= edge_radius_2.^2; 
% else    
%     edge_radius_2 = 907-2*(pp-1);  
%     edge_removal = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 > edge_radius_2.^2; 
% end

% figure,imshow(edge_removal)
TEMP_MAT = TEMP_MAT.*(edge_removal);

%% STEP 3: USE THE ACTUAL BINARIZED FEATURES TO REMOVE ANY UNWANTED FIBER FEATURES  THAT MAY BE PRESENT  
% LOAD THE H5 DATASET OF THE BINARIZED IMAGE VOLUME 
VOL = h5read(sprintf('%s%d%s%d%s%d%s','/scratch/halstead/r/ragyei/CRITICAL_REGIONS_CELL/',load_steps_idx(pp),'_files/XYZ_Recon_files_CONDENSED_N_',load_steps_idx(pp),'/Actual_Binary_Image_vol_1',load_steps_idx(pp),'_.h5'),'/CELL_DATA/FIBER_RECON');
% VOL_2 = VOL;
% VOL_2 = VOL_2.*uint8(edge_removal);
% 
% hh = 100;
% figure,imshow(logical(VOL(:,:,hh)))
% figure,imshow(logical(VOL_2(:,:,hh)))
% GET RID OF UNWANTED FIBER FEATURES USING THE BINARIZED IMAGE VOLUME AS A MASK
TEMP_MAT = TEMP_MAT.*double(VOL);


%% STEP 4: CREATE THE HDF OF THE MOST ACCURATE MASK
nom_file_prefix = 'OPTIMIZEDNET6666_';
[variable_out]  = MUTUAL_MAIN_HDF_CREATOR(TEMP_MAT,nom_file_prefix,load_steps_idx(pp),initial_directory_name);

% Remove the fibers on the outer edge and also remove the fiber fragments



%%
%     save(fullfile(initial_directory_name,sprintf('Idx_1_fin.mat')),'Idx_1_fin','-v7.3');
%     save(fullfile(initial_directory_name,sprintf('Idx_2_fin.mat')),'Idx_2_fin','-v7.3');
%     save(fullfile(initial_directory_name,sprintf('Idx_3_fin.mat')),'Idx_3_fin','-v7.3');

end

find(TEMP_MAT==304045);
numel(LIDX_var{142950})
unique(TEMP_MAT(LIDX_var{142950}))

hh = 142950

LIDX_var(304045)
[x,y,z] = ind2sub(size(TEMP_MAT),LIDX_var{304045});


% 
% for hh = 1:length(LIDX_var)
%     if all(TEMP_MAT(LIDX_var{hh}) == 0)
% 	   TEMP_MAT(LIDX_var{hh}) = hh;
%     else        
% %   break
%         unique_idx  =  unique(nonzeros(TEMP_MAT(LIDX_var{hh})));
%         tmp_idx     =  vertcat(unique_idx,hh);
%         LIDX_tmp    =  LIDX_var(tmp_idx);
% %         [~,~,z]     =  cellfun(@(x)ind2sub(size_length,x),LIDX_tmp,'uni',0);
%         z           =  cell2mat(cellfun(@(x)numel(unique(x)),z,'uni',0));    
%         if nnz(z == max(z)) == 1
%         [~,max_z] = max(z);
%         else
%         tmp_idx2 = find(z == max(z));    
%         [~,maxz2] = max(cell2mat(cellfun(@(x)numel(x),LIDX_tmp(z==max(z)),'uni',0))); 
%         max_z = tmp_idx2(maxz2);
%         end
%     
%     TEMP_MAT(vertcat(LIDX_tmp{:})) = 0;
%     LIDX_tmp(~ismember(1:numel(z),max_z)) = [];
%     TEMP_MAT(LIDX_tmp{1}) = hh;
%     
% %     % Update the variable LIDX_var
% %     [LIDX_var{z(~ismember(1:numel(z),max_z))}] = deal({});
%     end
% end
% 
% 
% 
% [x,y,z] = ind2sub(size_length,LIDX_var{252322});
% 
% [x,y,z] = ind2sub(size_length,LIDX{183100});



%     save(sprintf('%s%d%s','MULTI_VAR_',load_steps_idx(pp),'.mat'),'XYZ_ZXY_SM_SS','ZXY_XYZ_SM','empty_XYZ_ZXY',...
%                          'XYZ_YZX_SM_SS','YZX_XYZ_SM','empty_XYZ_YZX',...
%                          'ZXY_YZX_SM_SS','YZX_ZXY_SM','empty_ZXY_YZX',...
%                          'XYZ_ZXY_SM_size','XYZ_YZX_SM_size','ZXY_YZX_SM_size','-v7.3')

%     temp                      =     load(sprintf('%s%s',Center_Idx_prefix{sieve_indicator},'_xyz.mat'));
%     LIDX_XYZ                  =     temp.check_mat_center_GLOBAL_VEC_xyz_v2_UP;
%     temp                      =     load(sprintf('%s%s',CENT_COORD_cell{sieve_indicator},'_xyz.mat')); 
%     slice_XYZ                 =     cellfun(@(x)size(x,1),temp.FINAL_CORRECTED_CENTROID_xyz,'uni',0);


%     temp                      =     load(sprintf('%s%s',Center_Idx_prefix{sieve_indicator},'_zxy.mat'));
%     LIDX_ZXY                  =     temp.check_mat_center_GLOBAL_VEC_zxy_v2_UP;
%     temp                      =     load(sprintf('%s%s',CENT_COORD_cell{sieve_indicator},'_zxy.mat')); 
%     slice_ZXY                 =     cellfun(@(x)size(x,1),temp.FINAL_CORRECTED_CENTROID_zxy,'uni',0);

