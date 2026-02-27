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
% loop_vector_size_height    =    73;
sieve_indicator   =    1;
pp                =    1;

% % load the variable for the critical volume size
% load('loop_vector_size_height.mat')
% faulty_idx = [97 169 244];
% loop_vector_size_height = loop_vector_size_height(1:10,:);
% loop_vector_size_height(ismember(loop_vector_size_height(:,1),faulty_idx),:) = [];

loop_vector_size_height = [1 2]';

% loop_vector_size_height    =    loop_vector_size_height;
direc_prenom = '/scratch/halstead/r/ragyei/';
mfile_direc = '/scratch/halstead/r/ragyei/Mutual_Comparison_GITHUB_LINUX';

for pp = 1%:size(loop_vector_size_height,1)
    Center_Idx_prefix         =    {'Center_Idx_SIEVED','Center_Idx_UNSIEVED'};
%     initial_directory_name    =    sprintf('%s%d%s','/scratch/halstead/r/ragyei/CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files');  % Only make changes to this path 
    initial_directory_name    =    sprintf('%s%s%d%s',direc_prenom,'CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files');  % Only make changes to this path 

    
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
        cd(mfile_direc)        
        [REFINED_XYZ,REFINED_ZXY] = H5_Loader_Function(initial_directory_name,loop_vector_size_height(pp,1),sieve_indicator,REFINED_nom_cell,track_seq);
        size_length = size(REFINED_XYZ);
        
        cd(mfile_direc)
%         tic
        disp('begin EPAD_FUNCTION_CODE_MASTER_VERSION_MUT')
        [MAIN_COUNT_probe12,MAIN_COUNT_query12,SINGLE_CELL_12] = EPAD_FUNCTION_CODE_MASTER_VERSION_MUT(REFINED_XYZ,REFINED_ZXY);
%         toc

        clearvars   REFINED_XYZ       REFINED_ZXY 
        [Idx_1{ii},Idx_2{ii}] = MUTUAL_COMPARISON_CODE_ANALYZERS(MAIN_COUNT_probe12,MAIN_COUNT_query12,SINGLE_CELL_12);
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        IDX_PROBE = LIDX_XYZ;      IDX_QUERY = LIDX_ZXY; slice_PROBE = slice_XYZ;      slice_QUERY = slice_ZXY;     
       
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    elseif track_seq == 13
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    
        [REFINED_XYZ,REFINED_YZX] = H5_Loader_Function(initial_directory_name,loop_vector_size_height(pp,1),sieve_indicator,REFINED_nom_cell,track_seq);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cd(mfile_direc)
%         tic
        disp('begin EPAD_FUNCTION_CODE_MASTER_VERSION_MUT')
        [MAIN_COUNT_probe13,MAIN_COUNT_query13,SINGLE_CELL_13] = EPAD_FUNCTION_CODE_MASTER_VERSION_MUT(REFINED_XYZ,REFINED_YZX);
%         toc

        clearvars REFINED_XYZ      REFINED_YZX
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        [Idx_1{ii},Idx_3{ii}] = MUTUAL_COMPARISON_CODE_ANALYZERS(MAIN_COUNT_probe13,MAIN_COUNT_query13,SINGLE_CELL_13);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    elseif track_seq == 23
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        cd(mfile_direc)
        [REFINED_ZXY,REFINED_YZX] = H5_Loader_Function(initial_directory_name,loop_vector_size_height(pp,1),sieve_indicator,REFINED_nom_cell,track_seq);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        cd(mfile_direc)
%         tic
        disp('begin EPAD_FUNCTION_CODE_MASTER_VERSION_MUT')        
        [MAIN_COUNT_probe23,MAIN_COUNT_query23,SINGLE_CELL_23] = EPAD_FUNCTION_CODE_MASTER_VERSION_MUT(REFINED_ZXY,REFINED_YZX);
%         toc
      
        clearvars REFINED_ZXY     REFINED_YZX       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        [Idx_2{ii},Idx_3{ii}] = MUTUAL_COMPARISON_CODE_ANALYZERS(MAIN_COUNT_probe23,MAIN_COUNT_query23,SINGLE_CELL_23);
                        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    end
    end
        
Idx_1_fin   =   [Idx_1{:}];
Idx_2_fin   =   [Idx_2{:}];
Idx_3_fin   =   [Idx_3{:}];

% Clearvars to free up RAM
clearvars  MAIN_COUNT_probe12      MAIN_COUNT_probe13       MAIN_COUNT_probe23
clearvars  MAIN_COUNT_query12      MAIN_COUNT_query13       MAIN_COUNT_query23
clearvars  SINGLE_CELL_12          SINGLE_CELL_13           SINGLE_CELL_23
clearvars  Idx_1                   Idx_2                    Idx_3

%% Redefine the Linear Index of the fibers for the 3D volume:
LIDX       =   load(sprintf('%s%s%d%s%d%s',direc_prenom,'CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files/XYZ_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'/Center_Idx_UNSIEVED_xyz.mat'));
names      =   fieldnames(LIDX);
LIDX       =   LIDX.(names{1});
LIDX_var   =   LIDX(Idx_1_fin);

% clearvars LIDX 
% LIDX       =   load(sprintf('%s%d%s%d%s','/scratch/halstead/r/ragyei/CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files/XYZ_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'/Cent_Coord_UNSIEVED_xyz.mat'));
% names      =   fieldnames(LIDX);
% LIDX       =   LIDX.(names{1});
% Cent_var   =   LIDX(Idx_1_fin);

clearvars LIDX 
LIDX       =    load(sprintf('%s%s%d%s%d%s',direc_prenom,'CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files/ZXY_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'/Center_Idx_UNSIEVED_zxy.mat'));
names      =    fieldnames(LIDX);
LIDX       =    LIDX.(names{1});
LIDX_var   =   [LIDX_var     LIDX(Idx_2_fin)];

% clearvars LIDX 
% LIDX       =    load(sprintf('%s%d%s%d%s','/scratch/halstead/r/ragyei/CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files/ZXY_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'/Cent_Coord_UNSIEVED_zxy.mat'));
% names      =    fieldnames(LIDX);
% LIDX       =    LIDX.(names{1});
% Cent_var   =   [Cent_var     LIDX(Idx_2_fin)];


clearvars LIDX
LIDX       =   load(sprintf('%s%s%d%s%d%s',direc_prenom,'CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files/YZX_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'/Center_Idx_UNSIEVED_yzx.mat'));
names      =   fieldnames(LIDX);
LIDX       =   LIDX.(names{1});
LIDX_var   =   [LIDX_var     LIDX(Idx_3_fin)];

% clearvars LIDX
% LIDX       =   load(sprintf('%s%d%s%d%s','/scratch/halstead/r/ragyei/CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files/YZX_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'/Cent_Coord_UNSIEVED_yzx.mat'));
% names      =   fieldnames(LIDX);
% LIDX       =   LIDX.(names{1});
% Cent_var   =   [Cent_var     LIDX(Idx_3_fin)];


TEMP_MAT = zeros(size_length);
for hh = 1:length(LIDX_var)
    if all(TEMP_MAT(LIDX_var{hh}) == 0)
       if numel(LIDX_var{hh})>=100 
	   TEMP_MAT(LIDX_var{hh}) = hh;
       end
    else        
        unique_idx  =  unique(nonzeros(TEMP_MAT(LIDX_var{hh})));
        tmp_idx     =  vertcat(unique_idx,hh);
        LIDX_tmp    =  LIDX_var(tmp_idx);
        z           =  cell2mat(cellfun(@(x)numel(x),LIDX_tmp,'uni',0));    
        [~,max_z] = max(z);
    
    if tmp_idx(max_z) == hh    
    TEMP_MAT(vertcat(LIDX_tmp{:})) = 0;
    TEMP_MAT(LIDX_tmp{max_z}) = hh;  
    end
        
    end
end

% STEP 4: CREATE THE HDF OF THE MOST ACCURATE MASK
nom_file_prefix = 'OPTIMIZED_CRITICAL_VOL_PRIOR';
[variable_out]  = MUTUAL_MAIN_HDF_CREATOR(TEMP_MAT,nom_file_prefix,loop_vector_size_height(pp,1),initial_directory_name);


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

    edge_radius_2 = 0.97*centerX;
    edge_removal = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= edge_radius_2.^2; 

% figure,imshow(edge_removal)
TEMP_MAT = TEMP_MAT.*(edge_removal);

%% STEP 3: USE THE ACTUAL BINARIZED FEATURES TO REMOVE ANY UNWANTED FIBER FEATURES  THAT MAY BE PRESENT  
% LOAD THE H5 DATASET OF THE BINARIZED IMAGE VOLUME 
VOL = h5read(sprintf('%s%s%d%s%d%s%d%s',direc_prenom,'CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files/XYZ_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'/Actual_Binary_Image_vol_1',loop_vector_size_height(pp),'_.h5'),'/CELL_DATA/FIBER_RECON');
% VOL_2 = VOL;
% VOL_2 = VOL_2.*uint8(edge_removal);
% 
% hh = 100;
% figure,imshow(logical(VOL(:,:,hh)))
% figure,imshow(logical(VOL_2(:,:,hh)))
% GET RID OF UNWANTED FIBER FEATURES USING THE BINARIZED IMAGE VOLUME AS A MASK
TEMP_MAT = TEMP_MAT.*double(VOL);


%% STEP 4: CREATE THE HDF OF THE MOST ACCURATE MASK
nom_file_prefix = 'OPTIMIZED_CRITICAL_VOL_';
[variable_out]  = MUTUAL_MAIN_HDF_CREATOR(TEMP_MAT,nom_file_prefix,loop_vector_size_height(pp,1),initial_directory_name);

end


