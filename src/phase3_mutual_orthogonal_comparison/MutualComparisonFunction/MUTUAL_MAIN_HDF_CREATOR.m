function [varargout] = MUTUAL_MAIN_HDF_CREATOR(varargin)

TEMP_MAT            =    varargin{1};
nom_file_prefix     =    varargin{2};
new_idx             =    varargin{3};
save_path           =    varargin{4};


                  nom_bin            =   sprintf('%s%s%s',nom_file_prefix,new_idx,'_.h5');
                  nom_file           =   sprintf('%s%s%s',nom_file_prefix,new_idx,'_.xdmf');
                  delete (sprintf('%s',nom_bin))           
                  chunckSize         =   [200,200,50];           
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
           
           
                    h5att_name{1}      =    sprintf('%s%s%s',nom_file_prefix,'_xyz_',new_idx);
                    h5string{1}        =    sprintf('%s%s',nom_bin,':/CELL_DATA/FIBER_RECON');

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % ========================================================================
                    file_name          =   sprintf('%s',nom_file);
% % ========================================================================

                      z_depth      =   size(TEMP_MAT,3);
                      x_breadth    =   size(TEMP_MAT,2);
                      y_breadth    =   size(TEMP_MAT,1);
                      
                     
                XDMF_TEXT_GENERATOR(file_name,h5string,h5att_name,z_depth,y_breadth,x_breadth)   
                
%               newdir   =  sprintf('%s%d%s','G:\RECON_PROBE_FIX_NEW_v2\UNIQUE_FOLDER\RECON_',idx,'_files\XYZ_Recon_files');
                movefile(nom_bin,save_path)
                movefile(nom_file,save_path) 
                                              


stat_report = 'done';
varargout{    1} = stat_report;
