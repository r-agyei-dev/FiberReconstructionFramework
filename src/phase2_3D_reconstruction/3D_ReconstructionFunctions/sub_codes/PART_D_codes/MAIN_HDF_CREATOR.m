function [varargout] = MAIN_HDF_CREATOR(varargin)

TEMP_MAT            =    varargin{1};
nom_bin_prefix      =    varargin{2};
nom_file_prefix     =    varargin{3};
tomo_dataset_idx    =    varargin{4};
hdfsave             =    varargin{5};
load_variable       =    varargin{6};
save_path           =    varargin{7};

if hdfsave == 1    
     if load_variable == 1
                  nom_bin            =   sprintf('%s%s%s%s',nom_bin_prefix,'_xyz_',tomo_dataset_idx,'_.h5');
                  nom_file           =   sprintf('%s%s%s%s',nom_file_prefix,'_xyz_',tomo_dataset_idx,'_.xdmf');
                  delete (sprintf('%s',nom_bin))           
%                   chunckSize         =   [200,200, 50];      
                  chunckSize         =   round(0.25*size(TEMP_MAT));
                  h5create(sprintf('%s',nom_bin),'/CELL_DATA/FIBER_RECON',size(TEMP_MAT),'Datatype','double','ChunkSize',chunckSize,'Deflate',6 );                

                  [nx,ny,nz]         =   size(TEMP_MAT);
                  dz                 =   chunckSize(3);
                  zStart             =   1 : dz : nz;
                  zEnd               =   dz: dz : nz;
                  
                     if zStart(end) == nz
                        zStart  = zStart(1:end-1);
                     end
                     if zEnd(end) ~= nz
                        zEnd    = [zEnd,nz];
                     end
                 
                     for i = 1 : numel(zStart)
                           dz = zEnd(i) - zStart(i) + 1;
                             if dz > 0
                                start   =   [1 1 zStart(i)];                                                      % specifies the start of the dataset to start writing  
                                count   =   [nx ny dz];                                                 % numeric value specifying where in the data to start writing 
                                data    =   TEMP_MAT(:,:,zStart(i) :zEnd(i));                              % data to be written to the HDF5 file 
                                h5write(sprintf('%s',nom_bin),'/CELL_DATA/FIBER_RECON',data,start,count);
                             else
                                disp('dz == 0, cannot write data')
                             break
                             end
                     end
           
           
                    h5att_name{1}      =    sprintf('%s%s%s',nom_bin_prefix,'_xyz_',tomo_dataset_idx);
                    h5string{1}        =    sprintf('%s%s',nom_bin,':/CELL_DATA/FIBER_RECON');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ========================================================================
                    file_name          =   sprintf('%s',nom_file);
% % ========================================================================
                
       elseif  load_variable == 2  
                    
                  nom_bin            =   sprintf('%s%s%s%s',nom_bin_prefix,'_zxy_',tomo_dataset_idx,'_.h5');
                  nom_file           =   sprintf('%s%s%s%s',nom_file_prefix,'_zxy_',tomo_dataset_idx,'_.xdmf');
                  
%                   nom_bin            =   sprintf('%s%d%s','REFINED_FIBERS_zxy_',sieve_num,'.h5');
%                   nom_file           =   sprintf('%s%d%s','REFINED_FIBERS_zxy_',sieve_num,'.xdmf');

%                   nom_bin            =   'LUMP_REGION_zxy_new_app_UODATE.h5';
%                   nom_file           =   'LUMP_REGION_zxy_new_app_UODATE.xdmf';
                  delete (sprintf('%s',nom_bin))           
%                   chunckSize         =   [200,200, 50];    
                  chunckSize         =   round(0.25*size(TEMP_MAT));
                  h5create(sprintf('%s',nom_bin),'/CELL_DATA/FIBER_RECON',size(TEMP_MAT),'Datatype','double','ChunkSize',chunckSize,'Deflate',6 );                

                  [nx,ny,nz]         =   size(TEMP_MAT);
                  dz                 =   chunckSize(3);
                  zStart             =   1 : dz : nz;
                  zEnd               =   dz: dz : nz;
                  
                     if zStart(end) == nz
                        zStart  = zStart(1:end-1);
                     end
                     if zEnd(end) ~= nz
                        zEnd    = [zEnd,nz];
                     end
                 
                     for i = 1 : numel(zStart)
                           dz = zEnd(i) - zStart(i) + 1;
                             if dz > 0
                                start   =   [1 1 zStart(i)];                                                      % specifies the start of the dataset to start writing  
                                count   =   [nx ny dz];                                                 % numeric value specifying where in the data to start writing 
                                data    =   TEMP_MAT(:,:,zStart(i) :zEnd(i));                              % data to be written to the HDF5 file 
                                h5write(sprintf('%s',nom_bin),'/CELL_DATA/FIBER_RECON',data,start,count);
                             else
                                disp('dz == 0, cannot write data')
                             break
                             end
                     end
           
           
                    h5att_name{1}      =    sprintf('%s%s%s',nom_bin_prefix,'_zxy_',tomo_dataset_idx);
                    h5string{1}        =    sprintf('%s%s',nom_bin,':/CELL_DATA/FIBER_RECON');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ========================================================================
                    file_name          =   sprintf('%s',nom_file);
% % ========================================================================
    elseif  load_variable == 3     
        
                  nom_bin            =   sprintf('%s%s%s%s',nom_bin_prefix,'_yzx_',tomo_dataset_idx,'_.h5');
                  nom_file           =   sprintf('%s%s%s%s',nom_file_prefix,'_yzx_',tomo_dataset_idx,'_.xdmf');
                  
%                   nom_bin            =   'LUMP_REGION_yzx_new_app_UODATE.h5';
%                   nom_file           =   'LUMP_REGION_yzx_new_app_UODATE.xdmf';
                  delete (sprintf('%s',nom_bin))           
%                   chunckSize         =   [200,200, 50];    
                  chunckSize         =   round(0.25*size(TEMP_MAT));
                  h5create(sprintf('%s',nom_bin),'/CELL_DATA/FIBER_RECON',size(TEMP_MAT),'Datatype','double','ChunkSize',chunckSize,'Deflate',6 );                

                  [nx,ny,nz]         =   size(TEMP_MAT);
                  dz                 =   chunckSize(3);
                  zStart             =   1 : dz : nz;
                  zEnd               =   dz: dz : nz;
                  
                     if zStart(end) == nz
                        zStart  = zStart(1:end-1);
                     end
                     if zEnd(end) ~= nz
                        zEnd    = [zEnd,nz];
                     end
                 
                     for i = 1 : numel(zStart)
                           dz = zEnd(i) - zStart(i) + 1;
                             if dz > 0
                                start   =   [1 1 zStart(i)];                                                      % specifies the start of the dataset to start writing  
                                count   =   [nx ny dz];                                                 % numeric value specifying where in the data to start writing 
                                data    =   TEMP_MAT(:,:,zStart(i) :zEnd(i));                              % data to be written to the HDF5 file 
                                h5write(sprintf('%s',nom_bin),'/CELL_DATA/FIBER_RECON',data,start,count);
                             else
                                disp('dz == 0, cannot write data')
                             break
                             end
                     end
           
           
                    h5att_name{1}      =    sprintf('%s%s%s',nom_bin_prefix,'_yzx_',tomo_dataset_idx);
                    h5string{1}        =    sprintf('%s%s',nom_bin,':/CELL_DATA/FIBER_RECON');


                    
                    
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ========================================================================
                    file_name          =   sprintf('%s',nom_file);
% % ========================================================================   

                            
                
     end
                      z_depth      =   size(TEMP_MAT,3);
                      x_breadth    =   size(TEMP_MAT,2);
                      y_breadth    =   size(TEMP_MAT,1);
                      
                     
                XDMF_TEXT_GENERATOR(file_name,h5string,h5att_name,z_depth,y_breadth,x_breadth)   
                
%                 newdir   =  sprintf('%s%d%s','G:\RECON_PROBE_FIX_NEW_v2\UNIQUE_FOLDER\RECON_',idx,'_files\XYZ_Recon_files');
                movefile(nom_bin,save_path)
                movefile(nom_file,save_path) 
                                              
end

stat_report = 'done';
varargout{    1} = stat_report;
