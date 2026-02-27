function [varargout] = H5_Loader_Function(varargin)

    initial_directory_name    =    varargin{1};
    load_steps_idx_pp         =    varargin{2};
    sieve_indicator           =    varargin{3};
    REFINED_nom_cell          =    varargin{4};
    track_seq                 =    varargin{5};
       
    if track_seq == 12
     
%     E:\BACKUP_3_24_2021\specimen3_gfrp_0002\RECON_FILES\specimen3_gfrp_0002\XYZ_REC_FILES   
        
    disp('LOADING XYZ and ZXY')
    fileFolder                =    sprintf('%s%s%d',initial_directory_name,'/XYZ_Recon_files_CONDENSED_N_',load_steps_idx_pp);
    cd(fileFolder)
    tic
    disp('load curr')
    CURR_VOL               =     h5read(sprintf('%s%s%d%s',REFINED_nom_cell{sieve_indicator},'_xyz_',load_steps_idx_pp,'_.h5'),'/CELL_DATA/FIBER_RECON');
    toc

    fileFolder                =    sprintf('%s%s%d',initial_directory_name,'/ZXY_Recon_files_CONDENSED_N_',load_steps_idx_pp);
    cd(fileFolder)
    tic
    disp('load succ')
    SUCC_VOL               =     h5read(sprintf('%s%s%d%s',REFINED_nom_cell{sieve_indicator},'_zxy_',load_steps_idx_pp,'_.h5'),'/CELL_DATA/FIBER_RECON');
    toc
   
    elseif track_seq == 13
    disp('LOADING XYZ and YZX')       
    fileFolder                =    sprintf('%s%s%d',initial_directory_name,'/XYZ_Recon_files_CONDENSED_N_',load_steps_idx_pp);
    cd(fileFolder)
    tic
    disp('load curr')
    CURR_VOL              =     h5read(sprintf('%s%s%d%s',REFINED_nom_cell{sieve_indicator},'_xyz_',load_steps_idx_pp,'_.h5'),'/CELL_DATA/FIBER_RECON');
    toc

    fileFolder                =    sprintf('%s%s%d',initial_directory_name,'/YZX_Recon_files_CONDENSED_N_',load_steps_idx_pp);
    cd(fileFolder)
    tic
    disp('load succ')
    SUCC_VOL               =     h5read(sprintf('%s%s%d%s',REFINED_nom_cell{sieve_indicator},'_yzx_',load_steps_idx_pp,'_.h5'),'/CELL_DATA/FIBER_RECON');
    toc        
            
        
    elseif track_seq == 23
    disp('LOADING ZXY and YZX')   
    fileFolder                =    sprintf('%s%s%d',initial_directory_name,'/ZXY_Recon_files_CONDENSED_N_',load_steps_idx_pp);
    cd(fileFolder)
    tic
    disp('load curr')
    CURR_VOL               =     h5read(sprintf('%s%s%d%s',REFINED_nom_cell{sieve_indicator},'_zxy_',load_steps_idx_pp,'_.h5'),'/CELL_DATA/FIBER_RECON');
    toc

    fileFolder                =    sprintf('%s%s%d',initial_directory_name,'/YZX_Recon_files_CONDENSED_N_',load_steps_idx_pp);
    cd(fileFolder)
    tic
    disp('load succ')
    SUCC_VOL               =     h5read(sprintf('%s%s%d%s',REFINED_nom_cell{sieve_indicator},'_yzx_',load_steps_idx_pp,'_.h5'),'/CELL_DATA/FIBER_RECON');
    toc        
                 
    end
      
    
    varargout{             1}   =  CURR_VOL;
    varargout{         end+1}   =  SUCC_VOL;
    
end
