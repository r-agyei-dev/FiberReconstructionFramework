% This file evaluates the fiber length distribution of the sub_volume as well as computes the a metric for the Aij Tensor  

function [Fiber_length] = LENGTH_Aij_Tensor_ESTIMATOR_FUNCTION(centroid_cell_to_keep,size_length,lump,good_fibers)
% function [Fiber_length] = LENGTH_Aij_Tensor_ESTIMATOR_FUNCTION(Centroid_MID_cell,check_mat_center_GLOBAL_lump,size_length,lump)

%%
% ======================================================================================================
% ====>>> STEP 1: Initialize the various variables such as centroid lumping, vertical height extractions
% ======================================================================================================
if lump == 1
CENT_LUMP                                                                 =   [centroid_cell_to_keep{:}];
else
CENT_LUMP                                                                 =    centroid_cell_to_keep;
end


% FIB_PIXS                                                                  =   pix_to_keep;
% FIB_PIXS(cell2mat(cellfun(@(x)size(x,1),CENT_LUMP,'uni',0))<10)           =   [];
good_fibers(cell2mat(cellfun(@(x)size(x,1),CENT_LUMP,'uni',0))<10)        =   [];
CENT_LUMP(cell2mat(cellfun(@(x)size(x,1),CENT_LUMP,'uni',0))<10)          =   [];
% [~,~,Z_IDX]                                                               =   cellfun(@(x)ind2sub([size_length(1),size_length(2),size_length(3)],x),FIB_PIXS,'uni',0);
% Z_IDX                                                                     =   cellfun(@(x)unique(x),Z_IDX,'uni',0);

% ================================================
% ====>>> STEP 2: Compute the fiber length 
% ================================================
% Fiber_length   =    cell2mat(cellfun(@(x)sqrt(size(x,1)^2 + norm(x(1,:)-x(end,:))^2),CENT_LUMP,'uni',0));           % ====>>> Computes the out-of-plane angle 
Fiber_length   =    1.3*cell2mat(cellfun(@(x)norm(x(1,:)-x(end,:)),CENT_LUMP,'uni',0)); 


% ===========================================================================================================================
% ====>>> STEP 3: Merge the vertical height with the centroid pixel locations and compute the centroid for each fiber volume
% ===========================================================================================================================
% CENT_LUMP                                                                 =   cellfun(@(x,y)[x y],CENT_LUMP,Z_IDX,'uni',0);
% CENT_LUMP                                                                 =   centroid_cell_to_keep;
CENT_LUMP_MEDIAN                                                          =   cellfun(@(x)round(mean(x,1)),CENT_LUMP,'uni',0);
CENT_LUMP_MEDIAN_IDX                                                      =   cellfun(@(x)sub2ind([size_length(1),size_length(2),size_length(3)],x(2),x(1),x(3)),CENT_LUMP_MEDIAN,'uni',0);
cent_recon_ellipse                                                        =   [CENT_LUMP_MEDIAN_IDX{:}];

% ======================================================================================================================
% ====>>> STEP 4: Segment the volume into circular annulus volumes and extract the pixels locations within each annulus
% ======================================================================================================================
% Finding the centroid for each pixel location. Make concentric circles around the specimen 
% Finding the location of the fibers within the respective annulus 
clearvars imageSizeX imageSizeY columnsInImage rowsInImage circlePixels_1 circlePixels_2 n

z_depth                                                                   =   size_length(3);
imageSizeX                                                                =   size_length(2);
imageSizeY                                                                =   size_length(1);
[columnsInImage, rowsInImage]                                             =   meshgrid(1:imageSizeX, 1:imageSizeY);
centerX                                                                   =   size_length(2)/2;  
centerY                                                                   =   size_length(2)/2;   
radius                                                                    =   size_length(2)/2; 
n_divisions                                                               =   8;
temp_radius                                                               =   radius/sqrt(n_divisions);
circlePixels                                                              =   zeros(imageSizeX,imageSizeX,n_divisions);
equal_aea_pix_loc                                                         =   cell(1,n_divisions);
ring_idx                                                                  =   cell(1,n_divisions);

tic
for n = 1 : n_divisions
 
       temp_dark_region          =    zeros(imageSizeX);    
       circlePixels(:,:,n)       =    (rowsInImage - centerY).^2 + (columnsInImage - centerX).^2 <= (temp_radius*sqrt(n)).^2;
%        figure,imshow(circlePixels(:,:,n))

       if n == 1
           
       equal_aea_pix_loc{n}      =    find(repmat(circlePixels(:,:,n),[1,1,z_depth]));  
%        figure,imshow(circlePixels(:,:,n)) 
       
       elseif n > 1
    
       equal_aea_pix_loc{n}      =    find(repmat((circlePixels(:,:,n-1) - circlePixels(:,:,n))<0, [1,1,z_depth]));
  
       temp_dark_region((circlePixels(:,:,n-1) - circlePixels(:,:,n))<0)  =  1;
%        figure,imshow(temp_dark_region)

       end

       ring_idx{n}               =    find(ismember(cent_recon_ellipse,equal_aea_pix_loc{n}));
end
toc

% compile the centroids info according to their location in the annulus

       CENT_LUMP_RECOMPILED      =    cellfun(@(x)CENT_LUMP(x),ring_idx,'uni',0);

% =============================================================================
% ====>>> STEP 5: Compute the azimuth and elevations for the Aij Tensor terms 
% =============================================================================
      azimuth_cell               =    cell(1,n_divisions);
      elevation_cell             =    cell(1,n_divisions);


for nn    =    1   :   length(CENT_LUMP_RECOMPILED)

     rcell                       =    CENT_LUMP_RECOMPILED{nn};
     row_difference              =    cell2mat(cellfun(@(x)abs(x(1,2)-x(end,2)),rcell,'uni',0));
     column_difference           =    cell2mat(cellfun(@(x)abs(x(1,1)-x(end,1)),rcell,'uni',0)); 
     height_difference           =    cell2mat(cellfun(@(x)abs(x(1,end)-x(end,end)),rcell,'uni',0));    
     azimuth                     =    zeros(1,size(rcell,2));
     elevation                   =    zeros(1,size(rcell,2));

     for   tt   =  1  :  size(rcell,2)
     
            r      =    rcell{tt};
            
% qualify the values of the azimuth and elevations (bear in mind  rows is y and columns is x and the format is columns/row))

            if       (r(1,1) < r(end,1)) && (r(1,2)<r(end,2)) %=====>> Case 1(x1<x2 and y1<y2)  Quadrant 1 
    
               azimuth(tt)       =     90 - atand(column_difference(tt)/row_difference(tt));
               elevation(tt)     =     90 - atand(height_difference(tt)/sqrt(row_difference(tt)^2 + column_difference(tt)^2)); 
   
            elseif   (r(1,1) < r(end,1)) && (r(1,2)>r(end,2)) %=====>> Case 2(x1<x2 and y1>y2)  Quadrant 2 
    
               azimuth(tt)        =     90 + atand(column_difference(tt)/row_difference(tt));
               elevation(tt)      =     90 - atand(height_difference(tt)/sqrt(row_difference(tt)^2 + column_difference(tt)^2)); 

            elseif   (r(1,1) > r(end,1)) && (r(1,2)>r(end,2)) %=====>> Case 3(x1>x2 and y1>y2)  Quadrant 3
    
               azimuth(tt)         =     -(90 - atand(column_difference(tt)/row_difference(tt)));
               elevation(tt)       =     -(90 - atand(height_difference(tt)/sqrt(row_difference(tt)^2 + column_difference(tt)^2)));   
     
            elseif   (r(1,1) > r(end,1)) && (r(1,2)<r(end,2)) %=====>> Case 4(x1>x2 and y1<y2)  Quadrant 4
    
               azimuth(tt)         =     -(90 - atand(column_difference(tt)/row_difference(tt)));
               elevation(tt)       =     -(90 - atand(height_difference(tt)/sqrt(row_difference(tt)^2 + column_difference(tt)^2)));

            end

     end
 
               azimuth_cell{nn}    =      azimuth;
               elevation_cell{nn}  =      elevation;
end

% fill up the Aij matrix 

               mat_fill            =      [1 2 3  5 6 9];
             
               Aij_tensor_section                  =      cell(1,length(ring_idx));
               Aij                                 =      cell(1,length(ring_idx));
               Aij_tensor_section_mean             =      cell(1,length(ring_idx));
               
               for tt = 1:length(ring_idx)
    
               [Aij{tt}]                           =      Aij_Tensor_Estimator(azimuth_cell{tt},elevation_cell{tt});
               Aij_tensor_section{tt}              =      zeros(6,length(ring_idx{tt}));

                  for ll = 1:6
                    Aij_tensor_section{tt}(ll,:)   =       cell2mat(cellfun(@(x)x(mat_fill(ll)),Aij{tt},'uni',0));
                  end
                    
               Aij_tensor_section_mean{tt}         =       mean(Aij_tensor_section{tt},2);

               end

               Avg_Aij_across_sections             =       horzcat(Aij_tensor_section_mean{:});

% plot the Aij matrix
tt = 1;
eyy = cell2mat(cellfun(@(x)0.5*std(x(tt,:)),Aij_tensor_section,'uni',0));
yyy = 1:n_divisions;
figure, errorbar(yyy,Avg_Aij_across_sections(tt,:),eyy,'LineWidth',2)
ylim([0 1])

hold on 
tt = 2;
eyy = cell2mat(cellfun(@(x)0.5*std(x(tt,:)),Aij_tensor_section,'uni',0));
yyy = 1:n_divisions;
errorbar(yyy,Avg_Aij_across_sections(tt,:),eyy,'LineWidth',2)
ylim([0 1])

hold on 
tt = 3;
eyy = cell2mat(cellfun(@(x)0.5*std(x(tt,:)),Aij_tensor_section,'uni',0));
yyy = 1:n_divisions;
errorbar(yyy,Avg_Aij_across_sections(tt,:),eyy,'LineWidth',2)
ylim([0 1])

hold on 
tt = 4;
eyy = cell2mat(cellfun(@(x)0.5*std(x(tt,:)),Aij_tensor_section,'uni',0));
yyy = 1:n_divisions;
errorbar(yyy,Avg_Aij_across_sections(tt,:),eyy,'LineWidth',2)
ylim([0 1])

hold on 
tt = 5;
eyy = cell2mat(cellfun(@(x)0.5*std(x(tt,:)),Aij_tensor_section,'uni',0));
yyy = 1:n_divisions;
errorbar(yyy,Avg_Aij_across_sections(tt,:),eyy,'LineWidth',2)
ylim([0 1])

hold on 
tt = 6;
eyy = cell2mat(cellfun(@(x)0.5*std(x(tt,:)),Aij_tensor_section,'uni',0));
yyy = 1:n_divisions;
errorbar(yyy,Avg_Aij_across_sections(tt,:),eyy,'LineWidth',2)
ylim([0 1])

legend ('A11','A12','A13','A22','A23','A33')
title('plot of the orientation tensor')

xlabel('radial position')
ylabel('Average Orientation Tensor')

