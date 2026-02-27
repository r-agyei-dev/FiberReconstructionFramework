%%%%%%%%%%% THIS FUNCTION IS ONE THAT COMPUTES THE SIMILTANOUES CENTROID COMPARISON AND GIVES THE OUTPUT "SIEVE_PARAMETR_UNIQUE_CELL" AND "REMNANT_ELLIPSES"

% ===>>  THE "SIEVE_PARAMETER_UNIQUE_CELL IS THE VARIABLE THAT SIEVES INDEX OF MATCHING ELLIPSES ON SUBSEQUENT SLICES PERSISTING FORM THE CURRENT SLICE 
% ===>>  REMNANT_ELLIPSES IS THE VARAIBLE THAT HAS THE INDEX OF THE "UNMATCHED ELLIPSES" ON EACH SLICE

%%% BRIEF NOTES ON THE FUNCTION %%%
% THE FOLLOWING STEPS ARE CONSIDERED FOR THE FUNCTION BELOW AS:

% ===>> (1)  COMPUTE THE DISTANCES OF ALL THE ELLIPSES BETWEEN SLICES USING PDIST2. 
% ===>> (2)  REARRANGE THE DISTANCES SUCH THAT EACH ROW HAS THE FIRST 5 LEAST DISTANCES AND MAKE THE OTHERS ZERO. 
% ===>> (3)  FIND THE INDEX OF THOSE DISTANCES.
% ===>> (4)  EXTRACT THE CENTROID OF THE FIRST 5 LEAST FURTHER ELLIPSES. 
% ===>> (5)  MAKE A COMPARION OF THE VARIATION OF THEIR X AND Y CENTROID IN (%) FORM. 
% ===>> (6)  FIND THE INDEX OF THE LEAST PERCENTAGE CHANGE. 
% ===>> (7)  CREATE A MATRIX CALLED THE FINAL_SIEVE_PARAMETER WHICH WILL INCLUDE THE FOLLOWING 
% ===>> (7A) CURRENT ELLIPSE NUMBER AS COLUMN 1 
% ===>> (7B) PERCENTAGE CHANGE NUMBER AS COLUMN 2
% ===>> (7C) CADIDATE ELLIPSE NUMBER AS COLUMN 3
% ===>> (8)  COMPUTE A COMPETING CURENT ELLIPSES SCHEMATIC THAT COMPARES THE RELATIONSHIP BETWEEN CURRENT WITH THE SAME CANDIDATE ELLIPSES ON THE BASIS OF THE MINIMAL (%) CHANGE AND MAKE THE ROWS ALL ZERO. 

% function [sieve_parameter_unique_CELL,remnant_ellipses]    =   SIMIL_CENTROID_COMPARISON (numberOfImageFiles,SLICE_REGIONS_STATS_STRUCTURE,Centroid_cell)
function [sieve_parameter_unique_CELL,remnant_ellipses]    =   SIMIL_CENTROID_COMPARISON_UPDATED(numberOfImageFiles,SLICE_UPDATE,CENTROID_UPDATE)
%%
sieve_parameter_unique_CELL                                =   cell(1,numberOfImageFiles);
remnant_ellipses                                           =   cell(1,numberOfImageFiles);

for tt       =    1 : 10;%numberOfImageFiles-1
% clc;
% clear all variables within the loop 
clearvars dist_matrix dist_matrix_sort logical_idx cord_idx tmp_idx Cen_info_next Cen_info_current diff_array column_percentage_change
clearvars column_percent_change_cell row_percent_change_cell total_percent_change total_percent_change_cell tmp_idx_vector
clearvars tmp_idx_vector nt bint multiplet tmp_sort total_percent_change_columns min_val min_val_idx 
clearvars adel adel_vec sieve_parameter_unique sieve_parameter unique_regions

current     =  tt;
next        =  tt+1;

% find the lengths 
length_next                = length(SLICE_UPDATE{next});           % length of the ellipses on the next slice 
length_current             = length(SLICE_UPDATE{current});        % length of the ellipses on the current slice 

% ===>>> compute the distance of all ellipses between two slices 
dist_matrix                =     pdist2(CENTROID_UPDATE{current},CENTROID_UPDATE{next}); 

% ===>>> NOTE the rows is the number of regions in the previous and columns is the number of regions in the next slice 
dist_matrix_sort           =     dist_matrix';
ii = 1:length_next ;
dist_matrix_sort(ii,:)     =     sort(dist_matrix_sort(ii,:));

% ===>>> make non-zero the first five least distances and others zero
% dist_matrix_sort_keep      =     dist_matrix_sort';
dist_matrix_sort(6:end,:)  =     0;                                        
dist_matrix_sort           =     dist_matrix_sort';

logical_idx_tmp            =     cell(1,length_current);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ii = 1:length_current

   logical_idx_tmp{ii}     =     ismember(dist_matrix(ii,:),dist_matrix_sort(ii,:));  
   
   checker = find(logical_idx_tmp{ii});
   
   if numel(checker)>5
   % if its more than 5 then it suggests that the points are some sides on a circle so just pick one since those points are probably not the best  
   corrector = cell(1,5);
   for jj = 1:5
   corrector{jj} = find(ismember(dist_matrix(ii,:),dist_matrix_sort(ii,jj)));        
         if numel(corrector{jj})>1
                  corrector{jj}(2:end)=[];
         end   
   end
   corrector = [corrector{:}];
   
   logical_idx_tmp{ii}(checker(~ismember(checker,corrector))) = 0;
   end
   
   clearvars checker
   
end

% This ensures that the variables are unique terms
min_val_idx    =   cellfun(@(x,y)find(ismember(x,y)),num2cell(dist_matrix,2),num2cell(dist_matrix_sort(:,1),2),'uni',0);
min_val_idx    =   cell2mat(cellfun(@(x)x(1),min_val_idx,'uni',0));
min_val        =   dist_matrix_sort(:,1);

%% 

% Updated bit of the code 

%%
% create a matrix of final sieve parameters which will include the following ellipse number, ellipse comparison number and candidate ellipse 
% NOTE the last dimension must be analysed for competing sieve ellipses or competing current ellipses 
sieve_parameter                  =     [1:length_current;min_val';min_val_idx']';

% Compute that of the competing current ellipses 
[nt,~]                           =     histc(sieve_parameter(:,3),unique(sieve_parameter(:,3)));
unique_regions                   =     unique(sieve_parameter(:,3));
multiplet                        =     find(nt>1);

% Delete the "poorly competing regions on the basis of the min percentage change
% First compute index to delete
    adel                         =     cell(1,length(multiplet));
for ii = 1:length(multiplet)
    aa_temp                      =     sieve_parameter(:,3)== unique_regions(multiplet(ii));
    [~,aa_temp]                  =     min(sieve_parameter(aa_temp,2));
    aa_temp_check                =     find(sieve_parameter(:,3)== unique_regions(multiplet(ii)));
    aa_temp_check(aa_temp)       =     [];
    adel{ii}                     =     aa_temp_check;
end

% Concatenate the vectors 
adel                                  =     adel(~cellfun('isempty',adel));
adel_vec                              =     vertcat(adel{:});

% Create the sieve_parameter_unique_CELL
% if tt == 6 
%     break 
% else 
% end
sieve_parameter_unique                =     sieve_parameter;
sieve_parameter_unique([adel_vec],:)  =     0;
sieve_parameter_unique_CELL{tt}       =     sieve_parameter_unique;

% Create the remnant ellipses cell
remnant_ellipses{tt+1}                =     find(~ismember(1:length_next,sieve_parameter_unique_CELL{tt}(:,end)));

%%%%%%%%%%%%%%%%%%%%%%%%
sprintf('%s%d','END COMP LOOP',tt)
%%%%%%%%%%%%%%%%%%%%%%%%
end
remnant_ellipses{1}                   =     find(sieve_parameter_unique_CELL{1}(:,1)==0); 



%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% logical_idx                      =     vertcat(logical_idx_tmp{:});              % find the index of the first least distances 
% [cord_idx(:,1),cord_idx(:,2)]    =     find(logical_idx);                        % Note here that the first row is the current region and the second is the prospective regions idx 
% cord_idx                         =     sortrows(cord_idx);                         
% 
% % clearvars a
% % [a(:,1),a(:,2)] = hist(cord_idx(:,1),unique(cord_idx(:,1)));
% % find(a(:,1)~=5)
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % break the idx into cells 
% tmp_idx                          =     mat2cell(cord_idx(:,end)',1,repmat(5,1,length_current)); % we berak the chunks into the size of the current
% tmp_idx_columns                  =     cell2mat(cellfun(@(x) x',tmp_idx,'un',0));
% 
% ii = 1:length_current; 
% [Cen_info_next(:,1)]             =     deal(CENTROID_UPDATE{next}([tmp_idx{ii}],1));         % Gives the 5 array extracted centroid information for the next slice 
% [Cen_info_next(:,2)]             =     deal(CENTROID_UPDATE{next}([tmp_idx{ii}],end));       % Gives the 5 array extracted centroid information for the next slice 
% 
% % Use the kron replicate the indices of the first slice 
% Cen_info_current(:,1)            =     kron(CENTROID_UPDATE{current}(:,1),ones(5,1));
% Cen_info_current(:,2)            =     kron(CENTROID_UPDATE{current}(:,end),ones(5,1));
% 
% % Use the bsxfun to find the difference 
% diff_array                       =     abs(bsxfun(@minus,Cen_info_current,Cen_info_next));
% 
% % Use the bsxfun to find the percentage change 
% column_percent_change            =    (bsxfun(@rdivide,diff_array(:,1),Cen_info_next(:,1)))*100;
% row_percent_change               =    (bsxfun(@rdivide,diff_array(:,end),Cen_info_next(:,end)))*100;
% total_percent_change             =     column_percent_change + row_percent_change;
% 
% % Break the variable again into the respective chunks in a cell again
% % total_percent_change_cell = mat2cell(total_percent_change',1,repmat(5,1,length_current));
% total_percent_change_columns     =     reshape(total_percent_change,[5,length_current]);
% 
% % Find the idx and value of the least percent change 
% [min_val, min_val_idx]           =     min(total_percent_change_columns);
% 
% ii = 1:length_current;
% min_val_idx                      =     diag(tmp_idx_columns((min_val_idx(ii)),ii));
