function [varargout] = H5_Loader_FUNC(varargin)

    initial_directory_name    =    varargin{1};
    tomo_dataset_idx          =    varargin{2};
    specimen_name             =    varargin{3};
    REFINED_nom_cell          =    varargin{4};
    track_seq                 =    varargin{5};
     
    
    if track_seq == 12  
        
    disp('LOADING XYZ and ZXY')
    % fileFolder                =    sprintf('%s%s%s%s%s',initial_directory_name,'/RECON_FILES/',specimen_name,tomo_dataset_idx,'/XYZ_files');
    fileFolder                =    sprintf('%s%s',initial_directory_name,'/XYZ_files');
    cd(fileFolder)
    tic
    disp('load curr')
    CURR_VOL               =     h5read(sprintf('%s%s%s%s',REFINED_nom_cell{1},'_xyz_',tomo_dataset_idx,'_.h5'),'/CELL_DATA/FIBER_RECON');
    toc

    fileFolder                =    sprintf('%s%s',initial_directory_name,'/ZXY_files');    
    cd(fileFolder)
    tic
    disp('load succ')
    SUCC_VOL               =     h5read(sprintf('%s%s%s%s',REFINED_nom_cell{1},'_zxy_',tomo_dataset_idx,'_.h5'),'/CELL_DATA/FIBER_RECON');
    toc
   
    elseif track_seq == 13
    disp('LOADING XYZ and YZX')       
    fileFolder                =    sprintf('%s%s',initial_directory_name,'/XYZ_files');   
    cd(fileFolder)
    tic
    disp('load curr')
    CURR_VOL              =     h5read(sprintf('%s%s%s%s',REFINED_nom_cell{1},'_xyz_',tomo_dataset_idx,'_.h5'),'/CELL_DATA/FIBER_RECON');
    toc

%     fileFolder                =    sprintf('%s%s%d',initial_directory_name,'/YZX_Recon_files_CONDENSED_N_',load_steps_idx_pp);
    % fileFolder                =    sprintf('%s%s%s%s%s',initial_directory_name,'\RECON_FILES\',specimen_name,tomo_dataset_idx,'\YZX_REC_FILES');
    fileFolder                =    sprintf('%s%s',initial_directory_name,'/YZX_files');
    cd(fileFolder)
    tic
    disp('load succ')
    SUCC_VOL               =     h5read(sprintf('%s%s%s%s',REFINED_nom_cell{1},'_yzx_',tomo_dataset_idx,'_.h5'),'/CELL_DATA/FIBER_RECON');
    toc        
            
        
    elseif track_seq == 23
    disp('LOADING ZXY and YZX')   
    fileFolder                =    sprintf('%s%s',initial_directory_name,'/ZXY_files');
    cd(fileFolder)
    tic
    disp('load curr')
    CURR_VOL               =     h5read(sprintf('%s%s%s%s',REFINED_nom_cell{1},'_zxy_',tomo_dataset_idx,'_.h5'),'/CELL_DATA/FIBER_RECON');
    toc

    % fileFolder                =    sprintf('%s%s%s%s%s',initial_directory_name,'\RECON_FILES\',specimen_name,tomo_dataset_idx,'\YZX_REC_FILES');
    fileFolder                =    sprintf('%s%s',initial_directory_name,'/YZX_files');    
    cd(fileFolder)
    tic
    disp('load succ')
    SUCC_VOL               =     h5read(sprintf('%s%s%s%s',REFINED_nom_cell{1},'_yzx_',tomo_dataset_idx,'_.h5'),'/CELL_DATA/FIBER_RECON');
    toc        
                 
    end
      
    
    varargout{             1}   =  CURR_VOL;
    varargout{         end+1}   =  SUCC_VOL;
    
end
