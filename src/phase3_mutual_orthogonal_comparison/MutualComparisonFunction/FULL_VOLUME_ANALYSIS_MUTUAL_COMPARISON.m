                                                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                                                                            % This is a snippet of the code for the compariosn between 
                                                                            % two loading steps that uses Arithmetic operations only
                                                                            % BY: Ronald F Agyei.
                                                                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Description of the initial terms:
    %%%%% INITIAL_FIBER_LUMPING_SIEVE_FUNC_S3 %%%%%
    %%% INPUT VARIABLES %%%
    % VOL_PROBE is the 3D matrix volume of the fibers in the current configuration
    % VOL_QUERY is the 3D matrix volume of the fibers in the consecutive configuration
    % IDX_PROBE is the pixel indices of the fibers in the current configuration 
    % IDX_QUERY is the pixel indices of the fibers in the consecutive configuration
    % slice_PROBE is the number of conswecutive slices in the current configuration
    % slice_QUERY is the number of conswecutive slices in the consecutive configuration

    %%% OUTPUT VARIABLES %%%
    % MAIN_COUNT_empty is the variable that extracts empty sets through direct comparison
    % MAIN_COUNT_probe12 is the variable extracts one-to-many extracts between current and consecutive
    % MAIN_COUNT_query12 is the variable extracts one-to-many extracts between current and consecutive
    % MAIN_COUNT_single12 is the variable extracts one-to-many extracts between current and consecutive
    % FALSE_XYZ_ZXY_SM is the variable extracts one-to-many extracts between current and consecutive
    % FALSE_ZXY_XYZ_SM is the variable extracts one-to-many extracts between current and consecutive

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%% SHARING_VOLUME_FUNCTION_IMPROVED_S3 %%%%%
    %%% OUTPUT VARIABLES %%% 
    % MAIN_COUNT_probeTMP12 is the variable that matches one to many matches from probe to query
    % MAIN_COUNT_queryTMP12 is the variable that matches one to many matches from probe to query

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                            
    
    
    %   slice_instance = [1 2 3]; % xyz // zxy // yzx

    %%  Invoke the function that performs an initial sieve of the fibers (INITIAL_FIBER_LUMPING_SIEVE_FUNC) and eventually concatenates it into 
    % (a) one-to-many probe-query instance 
    % (b) many-to-one probe-query instance
    % (c) single instances 
    % (d) empty instances 

    %% Invoke the function that resolves: ==>>> (SHARING_VOLUME_FUNCTION_IMPROVED)
    % (a) distinct probe variables that share a common query
    % (b) distinct query variables that share a common probe 
    % (c) distinct large probe and query variables that share a common fiber ONE to MANY probe->query instance that share some vestige with a MANY to ONE probe->query 

    % if track_seq = 12, then we track xyz-->zxy;
    % if track_seq = 13, then we track xyz-->yzx;
    % if track_seq = 23, then we track zxy-->yzx;
                                                                            
                                                                            
                                                                            
MUTUAL_nom_cell   =   {'MUTUAL_RESULT_OPTIMIZED_SEIVEDn'};
REFINED_nom_cell  =   {'REFINED_FIBERS_SIEVED'};
TEMP_CRUDE_cell   =   'TEMP_CRUDE_FULL_';
CENT_COORD_cell   =   {'Cent_Coord_UNSIEVED'};
load_steps_idx    =   73;
sieve_indicator   =    1;
pp                =    1;

for pp = 1
    Center_Idx_prefix         =    {'Center_Idx_SIEVED','Center_Idx_UNSIEVED'};
    initial_directory_name    =     'C:\Users\armsb119pc4-user\Desktop\Three_D_code_UPDATE\UPDATED_SEG';
    fileFolder                =     sprintf('%s%s%d',initial_directory_name ,'\XYZ_Recon_files_',load_steps_idx(pp));
    cd(fileFolder)
    REFINED_XYZ               =     h5read(sprintf('%s%s%d%s',REFINED_nom_cell{sieve_indicator},'_xyz_3_',load_steps_idx(pp),'_.h5'),'/CELL_DATA/FIBER_RECON');
    temp                      =     load(sprintf('%s%s',Center_Idx_prefix{sieve_indicator},'_xyz.mat'));
    LIDX_XYZ                  =     temp.check_mat_center_GLOBAL_VEC_xyz_v2_UP;
    temp                      =     load(sprintf('%s%s',CENT_COORD_cell{sieve_indicator},'_xyz.mat')); 
    slice_XYZ                 =     cellfun(@(x)size(x,1),temp.FINAL_CORRECTED_CENTROID_xyz,'uni',0);

    fileFolder                =     sprintf('%s%s%d',initial_directory_name ,'\YZX_Recon_files_',load_steps_idx(pp));
    cd(fileFolder)
    REFINED_YZX               =     h5read(sprintf('%s%s%d%s',REFINED_nom_cell{sieve_indicator},'_yzx_3_',load_steps_idx(pp),'_.h5'),'/CELL_DATA/FIBER_RECON');
    temp                      =     load(sprintf('%s%s',Center_Idx_prefix{sieve_indicator},'_yzx.mat'));
    LIDX_YZX                  =     temp.check_mat_center_GLOBAL_VEC_yzx_v2_UP;
    temp                      =     load(sprintf('%s%s',CENT_COORD_cell{sieve_indicator},'_yzx.mat')); 
    slice_YZX                 =     cellfun(@(x)size(x,1),temp.FINAL_CORRECTED_CENTROID_yzx,'uni',0);

    fileFolder                =     sprintf('%s%s%d',initial_directory_name ,'\ZXY_Recon_files_',load_steps_idx(pp));
    cd(fileFolder)
    REFINED_ZXY               =     h5read(sprintf('%s%s%d%s',REFINED_nom_cell{sieve_indicator},'_zxy_3_',load_steps_idx(pp),'_.h5'),'/CELL_DATA/FIBER_RECON');
    temp                      =     load(sprintf('%s%s',Center_Idx_prefix{sieve_indicator},'_zxy.mat'));
    LIDX_ZXY                  =     temp.check_mat_center_GLOBAL_VEC_zxy_v2_UP;
    temp                      =     load(sprintf('%s%s',CENT_COORD_cell{sieve_indicator},'_zxy.mat')); 
    slice_ZXY                 =     cellfun(@(x)size(x,1),temp.FINAL_CORRECTED_CENTROID_zxy,'uni',0);

    %% initiate the while loop 
    clearvars MAIN_COUNT

    loop_track = [12 13 23];

    for ii = 1:3

        track_seq = loop_track(ii);

    if track_seq == 12
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%         VOL_PROBE = REFINED_XYZ;     VOL_QUERY = REFINED_ZXY;      
       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic
        [MAIN_COUNT_probe12,MAIN_COUNT_query12,SINGLE_CELL_12] = EPAD_FUNCTION_CODE_MASTER_VERSION_MUT(REFINED_XYZ,REFINED_ZXY);
        toc
        
        clearvars   REFINED_XYZ       REFINED_ZXY
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
 
%        IDX_PROBE = LIDX_XYZ;      IDX_QUERY = LIDX_ZXY; slice_PROBE = slice_XYZ;      slice_QUERY = slice_ZXY; 

       [MAIN_COUNT_probeTMP12,MAIN_COUNT_queryTMP12]    =   STEP_1B_SECONDARY_FIBER_LUMPING_FUNC(MAIN_COUNT_probe12,MAIN_COUNT_query12,SINGLE_CELL_12); 
       
       % Extract the multiple and the single terms from the MAIN_COUNT_probe and the MAIN_COUNT_query
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      function [varargout]  =  MUTUAL_COMPARISON_CODE_ANALYZERS(varargin)
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
%         XYZ_ZXY_SM_SS                   =        [MAIN_COUNT_probeTMP12,MAIN_COUNT_single12];
%         XYZ_ZXY_SM_size                 =        size(MAIN_COUNT_probeTMP12,2);
%         ZXY_XYZ_SM                      =        cellfun(@(x)[x(2,:); x(1,:); x(3,:)],MAIN_COUNT_queryTMP12,'uni',0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    elseif track_seq == 13
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        VOL_PROBE = REFINED_XYZ;     VOL_QUERY = REFINED_YZX;      

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic
        [MAIN_COUNT_probe13_TMP_1,MAIN_COUNT_query13_TMP_1,SINGLE_CELL_13] = EPAD_FUNCTION_CODE_MASTER_VERSION_MUT(VOL_PROBE,VOL_QUERY);
        toc
        
        cleavars REFINED_XYZ      REFINED_YZX
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        
%         IDX_PROBE = LIDX_XYZ;      IDX_QUERY = LIDX_YZX;    slice_PROBE = slice_XYZ;      slice_QUERY = slice_YZX; 

       [MAIN_COUNT_probeTMP13,MAIN_COUNT_queryTMP13]    =   STEP_1B_SECONDARY_FIBER_LUMPING_FUNC(MAIN_COUNT_probe13,MAIN_COUNT_query13,MAIN_COUNT_single13); 

       % Extract the multiple and the single terms from the MAIN_COUNT_probe and the MAIN_COUNT_query
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      function [varargout]  =  MUTUAL_COMPARISON_CODE_ANALYZERS(varargin)       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       
%         XYZ_YZX_SM_SS         =   [MAIN_COUNT_probeTMP13,MAIN_COUNT_single13];
%         XYZ_YZX_SM_size       =   size(MAIN_COUNT_probeTMP13,2);
%         YZX_XYZ_SM            =   cellfun(@(x)[x(2,:); x(1,:); x(3,:)],MAIN_COUNT_queryTMP13,'uni',0);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    elseif track_seq == 23
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%     
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        tic
        [MAIN_COUNT_probe23_TMP_1,MAIN_COUNT_query23_TMP_1,SINGLE_CELL_23] = EPAD_FUNCTION_CODE_MASTER_VERSION_MUT(REFINED_ZXY,REFINED_YZX);
        toc
        
        cleavars REFINED_ZXY     REFINED_YZX       
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
%         IDX_PROBE = LIDX_ZXY;      IDX_QUERY = LIDX_YZX;    slice_PROBE = slice_ZXY;      slice_QUERY = slice_YZX; 

        [MAIN_COUNT_probeTMP23,MAIN_COUNT_queryTMP23]    =   STEP_1B_SECONDARY_FIBER_LUMPING_FUNC(MAIN_COUNT_probe23,MAIN_COUNT_query23,MAIN_COUNT_single23); 

       % Extract the multiple and the single terms from the MAIN_COUNT_probe and the MAIN_COUNT_query
       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%      function [varargout]  =  MUTUAL_COMPARISON_CODE_ANALYZERS(varargin)       
       %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
          
%         ZXY_YZX_SM_SS        =   [MAIN_COUNT_probeTMP23,MAIN_COUNT_single23];
%         ZXY_YZX_SM_size      =   size(MAIN_COUNT_probeTMP23,2);
%         YZX_ZXY_SM           =   cellfun(@(x)[x(2,:); x(1,:); x(3,:)],MAIN_COUNT_queryTMP23,'uni',0);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

    end
    end

    save(sprintf('%s%d%s','MULTI_VAR_',load_steps_idx(pp),'.mat'),'XYZ_ZXY_SM_SS','ZXY_XYZ_SM','empty_XYZ_ZXY',...
                         'XYZ_YZX_SM_SS','YZX_XYZ_SM','empty_XYZ_YZX',...
                         'ZXY_YZX_SM_SS','YZX_ZXY_SM','empty_ZXY_YZX',...
                         'XYZ_ZXY_SM_size','XYZ_YZX_SM_size','ZXY_YZX_SM_size','-v7.3')

end
    