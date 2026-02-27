                                                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                            % This is a snippet of the code for the compariosn between 
                                                                            % two loading steps that uses Arithmetic operations only
                                                                            % BY: Ronald F Agyei.
                                                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                            
                                                                            % version update includes detection of void nucleation at the fiber tips
clear; close all; clc;  

MUTUAL_nom_cell   =   {'MUTUAL_RESULT_OPTIMIZED_SEIVEDn'};
REFINED_nom_cell  =   {'REFINED_FIBERS_SIEVED'};
TEMP_CRUDE_cell   =   'TEMP_CRUDE_FULL_';
CENT_COORD_cell   =   {'Cent_Coord_UNSIEVED'};
% loop_vector_size_height    =    73;
sieve_indicator   =    1;
pp                =    1;

% loop_vector_size_height = [1 14 16 18 20 22]';
loop_vector_size_height = [1]';

% direc_prenom = '/scratch/halstead/r/ragyei/';
% mfile_direc = '/scratch/halstead/r/ragyei/Mutual_Comparison_GITHUB_LINUX';

% direc_prenom = 'C:\Users\armsb119pc4-user\Desktop\';
% mfile_direc = 'C:\Users\armsb119pc4-user\Desktop\PhD_CODES\MUTUAL_COMPARSION_LINUX_VERS_1_2019';

direc_prenom = '/home/test/ragyei/';
mfile_direc = '/home/test/ragyei/CLUSTER_GITHUB_CODES/Mutual_Comparison_GITHUB_LINUX';



for pp = 1:size(loop_vector_size_height,1)
    Center_Idx_prefix         =    {'Center_Idx_SIEVED','Center_Idx_UNSIEVED'};
%     initial_directory_name    =    sprintf('%s%s%d%s',direc_prenom,'PhD_CODES\IMAGE_VOL\',loop_vector_size_height(pp,1),'_files');  % Only make changes to this path 
    initial_directory_name    =    sprintf('%s%s%d%s',direc_prenom,'CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files');  % Only make changes to this path 

    %% initiate the while loop 
    clearvars MAIN_COUNT

    loop_track  =   [12 13 23];
    Idx_1       =   cell(1,3);
    Idx_2       =   cell(1,3);
    Idx_3       =   cell(1,3);
    
    for ii = 1:3
    track_seq = loop_track(ii);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    if track_seq == 12
        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
%       VOL_PROBE = REFINED_XYZ;     VOL_QUERY = REFINED_ZXY;      
        cd(mfile_direc)        
        [REFINED_XYZ,REFINED_ZXY] = H5_Loader_Function(initial_directory_name,loop_vector_size_height(pp,1),sieve_indicator,REFINED_nom_cell,track_seq);
        size_length = size(REFINED_XYZ);
        
        cd(mfile_direc)
        disp('begin EPAD_FUNCTION_CODE_MASTER_VERSION_MUT')
        [MAIN_COUNT_probe12,MAIN_COUNT_query12,SINGLE_CELL_12] = EPAD_FUNCTION_CODE_MASTER_VERSION_MUT(REFINED_XYZ,REFINED_ZXY);

        clearvars   REFINED_XYZ       REFINED_ZXY 
        [Idx_1{ii},Idx_2{ii}] = MUTUAL_COMPARISON_CODE_ANALYZERS(MAIN_COUNT_probe12,MAIN_COUNT_query12,SINGLE_CELL_12);
        
%       IDX_PROBE = LIDX_XYZ;      IDX_QUERY = LIDX_ZXY; slice_PROBE = slice_XYZ;      slice_QUERY = slice_ZXY;           
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    elseif track_seq == 13

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    
        [REFINED_XYZ,REFINED_YZX] = H5_Loader_Function(initial_directory_name,loop_vector_size_height(pp,1),sieve_indicator,REFINED_nom_cell,track_seq);
        cd(mfile_direc)
        disp('begin EPAD_FUNCTION_CODE_MASTER_VERSION_MUT')
        [MAIN_COUNT_probe13,MAIN_COUNT_query13,SINGLE_CELL_13] = EPAD_FUNCTION_CODE_MASTER_VERSION_MUT(REFINED_XYZ,REFINED_YZX);

        clearvars REFINED_XYZ      REFINED_YZX
     
        [Idx_1{ii},Idx_3{ii}] = MUTUAL_COMPARISON_CODE_ANALYZERS(MAIN_COUNT_probe13,MAIN_COUNT_query13,SINGLE_CELL_13);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    elseif track_seq == 23

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        
        cd(mfile_direc)
        [REFINED_ZXY,REFINED_YZX] = H5_Loader_Function(initial_directory_name,loop_vector_size_height(pp,1),sieve_indicator,REFINED_nom_cell,track_seq);
        cd(mfile_direc)
        disp('begin EPAD_FUNCTION_CODE_MASTER_VERSION_MUT')        
        [MAIN_COUNT_probe23,MAIN_COUNT_query23,SINGLE_CELL_23] = EPAD_FUNCTION_CODE_MASTER_VERSION_MUT(REFINED_ZXY,REFINED_YZX);
      
        clearvars REFINED_ZXY     REFINED_YZX       

        [Idx_2{ii},Idx_3{ii}] = MUTUAL_COMPARISON_CODE_ANALYZERS(MAIN_COUNT_probe23,MAIN_COUNT_query23,SINGLE_CELL_23);
                        
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    end
    end
        
Idx_1_fin   =   [Idx_1{:}];
Idx_2_fin   =   [Idx_2{:}];
Idx_3_fin   =   [Idx_3{:}];

disp('Clearvars to free up RAM')
clearvars  MAIN_COUNT_probe12      MAIN_COUNT_probe13       MAIN_COUNT_probe23
clearvars  MAIN_COUNT_query12      MAIN_COUNT_query13       MAIN_COUNT_query23
clearvars  SINGLE_CELL_12          SINGLE_CELL_13           SINGLE_CELL_23
clearvars  Idx_1                   Idx_2                    Idx_3

%% Redefine the Linear Index of the fibers for the 3D volume:

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Redefine the Linear Index of the fibers for the 3D volume')
% LIDX       =   load(sprintf('%s%s%d%s%d%s',direc_prenom,'PhD_CODES\IMAGE_VOL\',loop_vector_size_height(pp,1),'_files\XYZ_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'\Center_Idx_UNSIEVED_xyz.mat'));
LIDX       =   load(sprintf('%s%s%d%s%d%s',direc_prenom,'CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files/XYZ_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'/Center_Idx_UNSIEVED_xyz.mat'));
names      =   fieldnames(LIDX);
LIDX       =   LIDX.(names{1});
tmp_labels =   unique(Idx_1_fin);
LIDX_var   =   LIDX(tmp_labels);

% Add a section for the end points of the centroids 
% CIDX       =   load(sprintf('%s%s%d%s%d%s',direc_prenom,'PhD_CODES\IMAGE_VOL\',loop_vector_size_height(pp,1),'_files\XYZ_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'\Cent_Coord_UNSIEVED_xyz.mat'));
CIDX       =   load(sprintf('%s%s%d%s%d%s',direc_prenom,'CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files/XYZ_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'/Cent_Coord_UNSIEVED_xyz.mat'));
names      =   fieldnames(CIDX);
% CIDX       =   CIDX.(names{1});
CIDX       =   cellfun(@(x)[x(:,2) x(:,1) x(:,3)],CIDX.(names{1}),'uni',0); % xyz is actually yxz is col/row/dep
CIDX_var   =   CIDX(tmp_labels);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% find(tmp_labels == 5312)
% LIDX_var(3879)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars LIDX 
% LIDX       =    load(sprintf('%s%s%d%s%d%s',direc_prenom,'PhD_CODES\IMAGE_VOL\',loop_vector_size_height(pp,1),'_files\ZXY_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'\Center_Idx_UNSIEVED_zxy.mat'));
LIDX       =    load(sprintf('%s%s%d%s%d%s',direc_prenom,'CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files/ZXY_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'/Center_Idx_UNSIEVED_zxy.mat'));
names      =    fieldnames(LIDX);
LIDX       =    LIDX.(names{1});
tmp_labels =    unique(Idx_2_fin);
LIDX_var   =    [LIDX_var     LIDX(tmp_labels)];

% Add a section for the end points of the centroids 
% CIDX       =   load(sprintf('%s%s%d%s%d%s',direc_prenom,'PhD_CODES\IMAGE_VOL\',loop_vector_size_height(pp,1),'_files\ZXY_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'\Cent_Coord_UNSIEVED_zxy.mat'));
CIDX       =   load(sprintf('%s%s%d%s%d%s',direc_prenom,'CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files/ZXY_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'/Cent_Coord_UNSIEVED_zxy.mat'));
names      =   fieldnames(CIDX);
% CIDX     =   CIDX.(names{1});
CIDX       =   cellfun(@(x)[x(:,1) x(:,3) x(:,2)],CIDX.(names{1}),'uni',0); % zxy is actually xzy is col/row/dep
CIDX_var   =   [CIDX_var      CIDX(tmp_labels)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LIDX_var(3879)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars LIDX
% LIDX       =   load(sprintf('%s%s%d%s%d%s',direc_prenom,'PhD_CODES\IMAGE_VOL\',loop_vector_size_height(pp,1),'_files\YZX_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'\Center_Idx_UNSIEVED_yzx.mat'));
LIDX       =   load(sprintf('%s%s%d%s%d%s',direc_prenom,'CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files/YZX_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'/Center_Idx_UNSIEVED_yzx.mat'));
names      =   fieldnames(LIDX);
LIDX       =   LIDX.(names{1});
tmp_labels =   unique(Idx_3_fin);
LIDX_var   =   [LIDX_var     LIDX(tmp_labels)];

% Add a section for the end points of the centroids 
% CIDX       =   load(sprintf('%s%s%d%s%d%s',direc_prenom,'PhD_CODES\IMAGE_VOL\',loop_vector_size_height(pp,1),'_files\YZX_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'\Cent_Coord_UNSIEVED_yzx.mat'));
CIDX       =   load(sprintf('%s%s%d%s%d%s',direc_prenom,'CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files/YZX_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'/Center_Idx_UNSIEVED_yzx.mat'));
names      =   fieldnames(CIDX);
% CIDX     =   CIDX.(names{1});
CIDX       =   cellfun(@(x)[x(:,3) x(:,2) x(:,1)],CIDX.(names{1}),'uni',0); % yzx is actually zyx is col/row/dep
CIDX_var   =   [CIDX_var      CIDX(tmp_labels)];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LIDX_var(3879)
%%

disp('Fill FULL_CRITICAL_VOL_UNLUMPED')
TEMP_MAT = zeros(size_length);
for hh = 1:length(LIDX_var)
	   TEMP_MAT(LIDX_var{hh}) = hh;
end

% STEP 4: CREATE THE HDF OF THE MOST ACCURATE MASK
disp('Create the HDF of the most accurate mask')
nom_file_prefix = 'FULL_CRITICAL_VOL_UNLUMPED';
[variable_out]  = MUTUAL_MAIN_HDF_CREATOR(TEMP_MAT,nom_file_prefix,loop_vector_size_height(pp,1),initial_directory_name);

%%
% Fill in fibers and replace based off the following 
% (a) Check if the there are no fibers and if the fiber length is above a threshold then fill 
% (b) If there competing overlapping fibers maintain the fibers with the longest fibers 

disp('Fill in fibers based off certain parameters')
tmp_idx = 0;
TEMP_MAT = zeros(size_length);
overlap_thresh = 0.2;

for hh = 1:length(LIDX_var)
    
    % check to see if the fiber is not a fragment and delete (may change this)
    if all(TEMP_MAT(LIDX_var{hh}) == 0)
%        if numel(LIDX_var{hh})>=100 
	   TEMP_MAT(LIDX_var{hh}) = hh;
%        end
    else
        % Isolate the overlapping indexes and extract the one with the largest fiber 
        unique_idx  =   unique(nonzeros(TEMP_MAT(LIDX_var{hh})));
        tmp_idx     =   vertcat(unique_idx,hh);
        LIDX_tmp    =   LIDX_var(tmp_idx);
        z           =   cell2mat(cellfun(@(x)numel(x),LIDX_tmp,'uni',0));    
        [~,max_z]   =   max(z);
        
        % Find the ratio of the overlap wrt the largest fiber
        LIDX_largest                   =     LIDX_tmp{max_z};
        LIDX_tmp(max_z)                =     [];
        tmp_idx_overlap                =     tmp_idx(~ismember(1:numel(tmp_idx),max_z));  % should have all indexes which are not the max
        numel_fib_pixels_tmp_overlap   =     cell2mat(cellfun(@(x)numel(x),LIDX_tmp,'uni',0));
        overlap_var                    =     cell2mat(cellfun(@(x)nnz(ismember(x,LIDX_largest))/numel(x),LIDX_tmp,'uni',0));
    
                % if the largest is the current fiber, remove all the smaller overlapping terms according to a threshold overlap
                if tmp_idx(max_z) == hh  
                LIDX_tmp = LIDX_tmp(overlap_var>overlap_thresh); 
                TEMP_MAT(vertcat(LIDX_tmp{:})) = 0;
                TEMP_MAT(LIDX_largest) = hh;  

                else
                    % if the largest is not the current fiber but needs to be filled in due to smaller overlap
                    if any(ismember(tmp_idx_overlap(overlap_var<overlap_thresh),hh))

                        if numel(tmp_idx_overlap)==1
                           TEMP_MAT(LIDX_var{hh})        =     hh; 
                        else 
                        % Apply the same sufficient overlap analogy to ONLY fibers smaller than the current fiber  

                                % Remove the current overlap idx in the current idxs
                                numel_hh                                               =     numel_fib_pixels_tmp_overlap(tmp_idx_overlap == hh);
                                numel_fib_pixels_tmp_overlap(tmp_idx_overlap == hh)    =     [];
                                tmp_idx_overlap(tmp_idx_overlap == hh)                 =     [];
                                LIDX_tmp_2                                             =     LIDX_var(tmp_idx_overlap);

                                % Look for the smaller ones
                                small_size_query = cell2mat(cellfun(@(x)numel(x),LIDX_tmp_2,'uni',0));

                                %                         numel_fib_pixels_tmp_overlap(small_size_query > numel(LIDX_var{hh}))    =     [];
                                %                         tmp_idx_overlap(small_size_query > numel(LIDX_var{hh}))                 =     [];
                                %                         numel_fib_pixels_tmp_overlap3     =    numel_fib_pixels_tmp_overlap(small_size_query > numel(LIDX_var{hh}));

                                tmp_idx_overlap3          =    tmp_idx_overlap(small_size_query > numel(LIDX_var{hh}));
                                LIDX_tmp_3                =     LIDX_var(tmp_idx_overlap3);
                                overlap_var_3             =    cell2mat(cellfun(@(x)nnz(ismember(x,LIDX_var{hh}))/numel(x),LIDX_tmp_3,'uni',0));

                                if all(overlap_var_3 < overlap_thresh)

                                    %                               numel_fib_pixels_tmp_overlap2     =    numel_fib_pixels_tmp_overlap(small_size_query <= numel(LIDX_var{hh}));

                                    tmp_idx_overlap2                                       =     tmp_idx_overlap(small_size_query <= numel(LIDX_var{hh}));
                                    LIDX_tmp_2                                             =     LIDX_var(tmp_idx_overlap2);

                                    if ~isempty(LIDX_tmp_2)
                                        overlap_var_2                                          =     cell2mat(cellfun(@(x)nnz(ismember(x,LIDX_var{hh}))/numel(x),LIDX_tmp_2,'uni',0));
                                        LIDX_tmp_2                                             =     LIDX_tmp_2(overlap_var_2>overlap_thresh);
                                        TEMP_MAT(vertcat(LIDX_tmp_2{:}))                       =     0;

                                        % fill in the new variable
                                        TEMP_MAT(LIDX_var{hh})                                 =     hh;
                                    end

                                end
            %             % find the numel of the current fibers 
            %             LIDX_tmp_2    =  LIDX_var(tmp_idx_overlap);
            %             z             =  cell2mat(cellfun(@(x)numel(x),LIDX_tmp_2,'uni',0));
            %             tmp_idx_overlap(z>numel(LIDX_var{hh})) = [];
            %             LIDX_tmp_2 = LIDX_var(tmp_idx_overlap);
            %             overlap_var_2 = cell2mat(cellfun(@(x)nnz(ismember(x,LIDX_var{hh}))/numel(x),LIDX_tmp_2,'uni',0));
            %             
            %             LIDX_tmp_2 = LIDX_tmp_2(overlap_var_2>0.2);
            %             TEMP_MAT(vertcat(LIDX_tmp_2{:})) = 0;
            %             % fill in the new variable
            %             TEMP_MAT(LIDX_var{hh}) = hh;
                        end
                    end
            %         
                end
        
    end
    %     
%     if any(tmp_idx == 3879)
%         break
%     end

end

% STEP 4: CREATE THE HDF OF THE MOST ACCURATE MASK
disp('Create the HDF of the most accurate mask')
nom_file_prefix = 'FINAL_OPTIMIZED_CRITICAL_VOL_UNLUMPED';
[variable_out]  = MUTUAL_MAIN_HDF_CREATOR(TEMP_MAT,nom_file_prefix,loop_vector_size_height(pp,1),initial_directory_name);


% % load in the file and check to see if the fibers are being loaded in correctly 
% tmp_1 = h5read(sprintf('%s%s',initial_directory_name,'/FINAL_OPTIMIZED_CRITICAL_VOL_UNLUMPED1_.h5'),'/CELL_DATA/FIBER_RECON');
% VizLayer(tmp_1,tmp_1)


% STEP 5: REMOVE ALL FIBER PROTRUDINGS
nom_file_prefix = 'FINAL_OPTIMIZED_CRITICAL_VOL_UNLUMPED_NO_PROTRUDE';
PROTRUDE_CORRECTOR = h5read(sprintf('%s%s%d%s%d%s%d%s',direc_prenom,'CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files/XYZ_Recon_files_CONDENSED_N_',loop_vector_size_height(pp,1),'/Actual_Binary_Image_vol_1',loop_vector_size_height(pp,1),'_.h5')...
                            ,'/CELL_DATA/FIBER_RECON');
TEMP_MAT(PROTRUDE_CORRECTOR~=1) = 0;
[variable_out]  = MUTUAL_MAIN_HDF_CREATOR(TEMP_MAT,nom_file_prefix,loop_vector_size_height(pp,1),initial_directory_name);


%%

numel_idx = unique(nonzeros(TEMP_MAT));
DataOut1  = zeros(length(numel_idx),9);
tic 
    for i=1:length(numel_idx)                                                                                      %    for each fiber       

        % Use the variable  
        [X,Y,Z]              =     ind2sub(size_length,LIDX_var{numel_idx(i)});          %    the subscripts of the given indices  
        ZminIn               =     find(Z==min(Z));                                                             %    the indices of the minimum z
        ZmaxIn               =     find(Z==max(Z));                                                             %    the indices of the maximum z       
        startcord            =     [mean(X(ZminIn)) mean(Y(ZminIn)) mean(Z(ZminIn))];                           %    the x y and z cord of 'centroid'
        endcord              =     [mean(X(ZmaxIn)) mean(Y(ZmaxIn)) mean(Z(ZmaxIn))];                           %    the x y and z cord of 'centroid'
        
        % Eliminate disks:
%         if or( or( startcord(1)==endcord(1) , startcord(2)==endcord(2) )==1 , startcord(3)==endcord(3))==1 % if any of the enpoints equal each other in a direction
        if or( or( min(X)==max(X) , min(Y)==max(Y) )==1 , min(Z)==max(Z))==1 % if any of the enpoints equal each other in a direction
            
            DataOut1(i,:)=zeros(1,9); % do nothing, plug in all zeros in the output variable
        else % fiber is longer than one slice
            DataOut1(i,1:3)=startcord; % plug in coordinates for the starting endpoint
            DataOut1(i,4:6)=endcord; % plug in the coordinates for the ending endpoint
            % Angles:
            L=sqrt( (endcord(1)-startcord(1))^2 + (endcord(2)-startcord(2))^2 + (endcord(3)-startcord(3))^2 ); %the length of the fiber
            phi = acosd( (endcord(1)-startcord(1)) / L); % Phi, range is between 0 and 180
            theta = sign( endcord(2)-startcord(2) ) * acosd( (endcord(3)-startcord(3)) / L); % theta, range is between 0 and 90
            % plug in angles to the output file
            DataOut1(i,7) = phi;
            DataOut1(i,8) = theta ;
            DataOut1(i,9) = L;
        end
    end
toc 


    for i=1:length(numel_idx)                                                                                      %    for each fiber       

        % Use the variable  
        [X,Y,Z]                   =     ind2sub(size_length,LIDX_var{numel_idx(i)});          %    the subscripts of the given indices  
        centroid_variable(:,i)    =     mean([X,Y,Z]);
        
    end



end

%%

% save('centroid_variable.mat','centroid_variable','-v7.3')
% save('DataOut1.mat','DataOut1','-v7.3')

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % STEPS for ensuring accurate 3D fiber reconstructions 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % %% STEP 1: BINARIZE IMAGE VOLUME AND PERFORM SLICE-WISE REMOVAL OF ALL THE SMALL FEATURES ALONG ALL THE DIMENSIONS 
% % % BINARIZE 
% % TEMP_BIN = TEMP_MAT ~= 0;
% % 
% % % SMALL FEATURES REMOVAL ALONG THE XYZ DXN
% % for ii = 1:size(TEMP_BIN,3)    
% % TEMP_BIN(:,:,ii)  = bwareaopen(TEMP_BIN(:,:,ii),5);      
% % end
% % 
% % % SMALL FEATURES REMOVAL ALONG THE ZXY DXN
% % TEMP_BIN                =   permute(TEMP_BIN,[3 1 2]);              % (z x y)
% % for ii = 1:size(TEMP_BIN,3)    
% % TEMP_BIN(:,:,ii)  = bwareaopen(TEMP_BIN(:,:,ii),5);      
% % end
% % TEMP_BIN              =   permute(TEMP_BIN,[2 3 1]);                % (z x y) ==>> xyz
% % % 
% % % SMALL FEATURES REMOVAL ALONG THE YZX DXN
% % TEMP_BIN                =   permute(TEMP_BIN,[2 3 1]);              % (y z x)
% % for ii = 1:size(TEMP_BIN,3)    
% % TEMP_BIN(:,:,ii)  = bwareaopen(TEMP_BIN(:,:,ii),5);      
% % end
% % TEMP_BIN               =   permute(TEMP_BIN,[3 1 2]);               % (y z x) ==>> xyz
% % %
% % 
% % % REMOVE THE 3D USING THE SLICE REMOVAL MASK
% % TEMP_MAT = TEMP_MAT.*TEMP_BIN;
% % 
% % %% STEP 2: REMOVE ANY FIBER THAT MAY BE WITHIN THE BAND OF FALSE RECONSTRUCTION AROUND THE SPECIMEN 
% % 
% % imageSizeX = size(TEMP_MAT,2);  imageSizeY = size(TEMP_MAT,1);
% % [columnsInImage,rowsInImage] = meshgrid(1:imageSizeX, 1:imageSizeY);
% % % Next create the circle in the image.
% % centerX = 0.5*imageSizeX;  centerY = 0.5*imageSizeY;
% % 
% %     edge_radius_2 = 0.97*centerX;
% %     edge_removal = (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= edge_radius_2.^2; 
% % 
% % % figure,imshow(edge_removal)
% % TEMP_MAT = TEMP_MAT.*(edge_removal);
% % 
% % %% STEP 3: USE THE ACTUAL BINARIZED FEATURES TO REMOVE ANY UNWANTED FIBER FEATURES  THAT MAY BE PRESENT  
% % % LOAD THE H5 DATASET OF THE BINARIZED IMAGE VOLUME 
% % VOL = h5read(sprintf('%s%s%d%s%d%s%d%s',direc_prenom,'CRITICAL_REGIONS_CELL/',loop_vector_size_height(pp,1),'_files/XYZ_Recon_files_CONDENSED_N_',loop_vector_size_height(pp),'/Actual_Binary_Image_vol_1',loop_vector_size_height(pp),'_.h5'),'/CELL_DATA/FIBER_RECON');
% % % VOL_2 = VOL;
% % % VOL_2 = VOL_2.*uint8(edge_removal);
% % % 
% % % hh = 100;
% % % figure,imshow(logical(VOL(:,:,hh)))
% % % figure,imshow(logical(VOL_2(:,:,hh)))
% % % GET RID OF UNWANTED FIBER FEATURES USING THE BINARIZED IMAGE VOLUME AS A MASK
% % TEMP_MAT = TEMP_MAT.*double(VOL);
% % 
% % 
% % %% STEP 4: CREATE THE HDF OF THE MOST ACCURATE MASK
% % nom_file_prefix = 'OPTIMIZED_CRITICAL_VOL_';
% % [variable_out]  = MUTUAL_MAIN_HDF_CREATOR(TEMP_MAT,nom_file_prefix,loop_vector_size_height(pp,1),initial_directory_name);
