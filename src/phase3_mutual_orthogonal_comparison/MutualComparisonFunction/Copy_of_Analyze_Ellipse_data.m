close all
clear all
clc

% ed=csvread('S1_Results_layer1.csv');
% 
% % image_name='S3_20x_circular.tif';
% image_name='S1_Full_Etched.tif';
% 
% image=imread(image_name);


%%
load('DataOut1.mat')
load('centroid_variable.mat')
centroid_variable = centroid_variable';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %S1 layer 1:
% center=[5097 5145]; %S1 center is manually inputted [y x]

%S1 Layer 2
% center=[10421 10860]; %center is manually inputted [y x]
center = [946 946]; %center is manually inputted [y x]
radius_length = 1.3*sqrt((centroid_variable(:,1)-946).^2  + (centroid_variable(:,2)-946).^2 );
% %S3 layer 1:
% center=[9900 10100]; %S3 center is manually inputted [y x]

% % %For S3 layer 2
% center=[11700 11494]; %center is manually inputted [y x]


% % %For S4 layer 2
% center=[11300 9860]; %center is manually inputted [y x]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% diam=2.42; %S3 and S4 diameter of interest in mm, from tensile experiments
diam=2.46; %S3 and S4 diameter of interest in mm, from tensile experiments
% diam=2.48; %S1 diameter of interest in mm, from tensile experiments


% ellipseDATA2{i} = [xc{i} yc{i} major{i} minor{i} phi{i} theta{i} radiusn{i}'];

xc= centroid_variable(:,1); %ed(:,1);
yc= centroid_variable(:,2); %ed(:,2);
% major=ed(:,3);
% minor=ed(:,4);
phi=DataOut1(:,7); %ed(:,5);
theta=DataOut1(:,8); %ed(:,6);
% r=ed(:,7);
r = radius_length;
A11=(sind(theta).^2).*(cosd(phi).^2);
A22=(sind(theta).^2).*(sind(phi).^2);
A33=(cosd(theta).^2);
A12=(sind(theta).^2).*(cosd(phi)).*(sind(phi));

sections=10;
% rm=max(r)+4;
% rm=diam/2*1000*1000/323.09;
rm = 1.3*946;

% % Figured out each radius of the total sections
for j=1:sections
    rs(j)=sqrt(j/sections*rm^2);
    centers(j,:)=[center(2) center(1)];
end

image = ones(1892);
figure
imshow(image)
hold on
viscircles(centers,rs)

for j=1:sections;
  k=1;
  
  if j==1
    for i=1:length(r);
        if r(i)<rs(1)
            A11_sections{j}(k)=A11(i); 
            A22_sections{j}(k)=A22(i); 
            A33_sections{j}(k)=A33(i); 
            A12_sections{j}(k)=A12(i);
            k=k+1;
        end    
    end
    center_radius(j)=rs(1)/2;
    mean_A11(j)=mean(A11_sections{j});
    std_A11(j)=std(A11_sections{j})/2;
    
    mean_A22(j)=mean(A22_sections{j});
    std_A22(j)=std(A22_sections{j})/2;
    
    mean_A33(j)=mean(A33_sections{j});
    std_A33(j)=std(A33_sections{j})/2;
    
    mean_A12(j)=mean(A12_sections{j});
    std_A12(j)=std(A12_sections{j})/2;
    
    length(A11_sections{j})
  elseif j>1
      for i=1:length(r);
        if r(i)>=rs(j-1) && r(i)<rs(j)
            A11_sections{j}(k)=A11(i); 
            A22_sections{j}(k)=A22(i); 
            A33_sections{j}(k)=A33(i); 
            A12_sections{j}(k)=A12(i); 
            k=k+1;
        end    
    end
    center_radius(j)=rs(j)-(rs(j)-rs(j-1))/2;
    mean_A11(j)=mean(A11_sections{j});
    std_A11(j)=std(A11_sections{j})/2;
    
    mean_A22(j)=mean(A22_sections{j});
    std_A22(j)=std(A22_sections{j})/2;
    
    mean_A33(j)=mean(A33_sections{j});
    std_A33(j)=std(A33_sections{j})/2;
    
    mean_A12(j)=mean(A12_sections{j});
    std_A12(j)=std(A12_sections{j})/2;
    length(A11_sections{j})
  end
  
end

figure
% plot(center_radius.*323.09/1000/1000,mean_A11, '*-')
% errorbar(center_radius.*323.09/1000/1000, mean_A11, std_A11, '*-')
errorbar(center_radius,mean_A11, std_A11, '*-')
ylim([0 1])
xlabel('Radius (mm)')
ylabel('A11')
title('S3 A11 distribution')

figure
% plot(center_radius.*323.09/1000/1000,mean_A11, '*-')
% errorbar(center_radius.*323.09/1000/1000, mean_A22, std_A22, '*-')
errorbar(center_radius,mean_A22, std_A22, '*-')
ylim([0 1])
xlabel('Radius (mm)')
ylabel('A22')
title('S3 A22 distribution')

figure
% plot(center_radius.*323.09/1000/1000,mean_A11, '*-')
% errorbar(center_radius.*323.09/1000/1000, mean_A33, std_A33, '*-')
errorbar(center_radius,mean_A33, std_A33, '*-')
ylim([0 1])
xlabel('Radius (mm)')
ylabel('A33')
title('S3 A33 distribution')

figure
% plot(center_radius.*323.09/1000/1000,mean_A11, '*-')
% errorbar(center_radius.*323.09/1000/1000, mean_A12, std_A12, '*-')
errorbar(center_radius,mean_A33, std_A33, '*-')
ylim([0 1])
xlabel('Radius (mm)')
ylabel('A12')
title('S3 A12 distribution')


% figure
% errorbar(center_radius.*323.09/1000/1000, mean_A11, std_A11, '*-')
% hold on
% errorbar(center_radius_S3.*323.09/1000/1000, mean_A11_S3, std_A11, '*-')
% ylim([0 1])
% xlabel('Radius (mm)')
% ylabel('A11')
% title('S1 and S3 A11 distribution')
% legend('S1', 'S3')

% All together now!

image_name='Aij Tensor Distribution';

figure
% errorbar(center_radius.*323.09/1000/1000, mean_A11, std_A11, '*-')
errorbar(center_radius, mean_A11, std_A11, '*-')
hold on
% errorbar(center_radius.*323.09/1000/1000, mean_A22, std_A11, '*-')
errorbar(center_radius, mean_A22, std_A11, '*-')
% errorbar(center_radius.*323.09/1000/1000, mean_A33, std_A11, '*-')
errorbar(center_radius, mean_A33, std_A11, '*-')
% ylim([0 1])
% xlim([0 1.24])
xlabel('Radius (mm)')
ylabel('A_{ij}')
title([image_name, ' A_{ij} distribution'])
legend('A11', 'A22', 'A33')
        