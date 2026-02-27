%% FUNCTION DESCRIPTION:
% This function analyzes the bi-directional stitching of the fibers. 
% This takes into acccount:

% (a) Mitigation of erroneous reconstruction of fiber segments with collinear longitudinal axis 
% separated by gaps due to subliminal intensity below detection of supervised iteration. 
% (b) Mitigation of collinear stacked fibers after the fiber separating algorithm. 
% (c) Mitigation of erroneous fiber separation due to in-plane over-segmentation of elliptical regions. 
% (d) Mitigation of over-segmented reconstructed fibers due to rogue regions.

% Code has been updated to ensure the use of the condensed format 

% The basic sieving approach includes an inplane distance threshold and an out of plane angle.
function [varargout] = THREE_D_STITCHER_OF_SEG_FIBERS_FUNCTION(varargin)

% Input variables allocation:
avail_terms                          =      varargin{1};                    % Index of the labelled fibers 
Linear_Index_center_GLOBAL_LUMP      =      varargin{2};                    % Linear_Idx of the labelled fibers 
size_length                          =      varargin{3};                    % size of the 3D volume 
minRadius_idxval_LUMP_cat            =      varargin{4};                    % Average diameter and numel of constituent ellipses for the labelled fibers 
Centroid_MID_cell_LUMP_ORIG          =      varargin{5};

clearvars    Centers               radius                            imageSize_row                         imageSize_col               columnsInImage         rowsInImage
clearvars    circlePixels_top      circlePixels_top_pt               Centroid_MID_cell_LUMP_UPDATE         grad_out_LUMP_UPDATE        temp                   fib_diam 
clearvars    likely_term_var       UPDATE_CHECK_MAT_GLOBAL_LUMP      z                            

%% STEP (I): Extract the beginning and terminal slices for the respective fibers
% Extract the pixel locations of the linear indices of the labelled fibers 
[~,~,ze]                     =     cellfun(@(x)ind2sub(size_length,x),Linear_Index_center_GLOBAL_LUMP,'uni',0);   % Chnge Here 

clearvars xe ye 

% Set TEMPORARY VALUES to account for the broadcast variables
minThresh    =  5;                      %   threshold for distance of centroids in stitching

% Run the condensed format to extract the respective fiber indices 
% Extracted information includes
% (a) Fiber end points/ fiber orientation of angle(in plane and out of plane)
% (b) Extract the code for the initial lumpoed fibers 
% (c) Review the interpolation and extrapolation schemes for the stitched fibers 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Begin analyzing
     
    disp('Condensing 3D volume data: gathering fiber end points and computing angles');  % display for the user
    fibers             =   avail_terms;
    DataOut1           =   zeros(length(fibers), 8);    %  initialize output data variable  
    
    for i=1:length(fibers)                                                                                      %    for each fiber       

        % Compute the startcord and endcord 
        [X,Y,Z]              =     ind2sub(size_length,Linear_Index_center_GLOBAL_LUMP{fibers(i)});             %    the subscripts of the given indices  
        ZminIn               =     find(Z==min(Z));                                                             %    the indices of the minimum z
        ZmaxIn               =     find(Z==max(Z));                                                             %    the indices of the maximum z       
        startcord            =     [mean(X(ZminIn)) mean(Y(ZminIn)) mean(Z(ZminIn))];                           %    the x y and z cord of 'centroid'
        endcord              =     [mean(X(ZmaxIn)) mean(Y(ZmaxIn)) mean(Z(ZmaxIn))];                           %    the x y and z cord of 'centroid'
        
        % Eliminate disks:
        % Create a variable of codnesed data that has the metadata of ellipses 
        if or( or( startcord(1)==endcord(1) , startcord(2)==endcord(2) )==1 , startcord(3)==endcord(3))==1 % if any of the enpoints equal each other in a direction
            DataOut1(i,:)=zeros(1,8); % do nothing, plug in all zeros in the output variable
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
        end
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Perform fiber stitching based of the endpoint proximity and the minimization of Euclidean angles
% Exclude fibers that have been seperated using the Fiber Seperating Algorithm (DO NOT RESTITCH)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Determining discontinuous fibers
%   Please ensure that always there is a unique match and no other cases of ambiguity  
%   issue_case = zeros(length(fibers),1);
    fib_stitch_thresh = 11;
    fib_stitch_thresh_length = 13;

      disp('Stitching discontinuous fibers in condensed data')
    j=1; % counter for stitched fibers
    for i=1:length(fibers)
%         sprintf('%s%d','fiber sticth ',i)
        diff = abs(bsxfun(@minus,DataOut1(i,4:6),DataOut1(:,1:3)));   %  difference of the end cord to start cord
        MinDifInX = find(diff(:,1)<minThresh);           %  find indices of fibers that are close in x
        MinDifInY = find(diff(:,2)<minThresh);           %  find indices of fibers that are close in y
        MinDifInZ = find(diff(:,3)<minThresh);           %  find indices of fibers that are close in z
        MinXY=intersect(MinDifInX, MinDifInY);           %  close in the x AND y
        MinXYZ=intersect(MinXY,MinDifInZ);               %  indices of fibers that are close in x AND y AND z (DO NOT CONFUSE INDICES WITH LABELS)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
        if DataOut1(i,1:8)==zeros(1,8)  %  if this was a disk that was deleted                                                              ======>>>> UNIQUE VERIFIED
            continue %skip it
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
        elseif or(isempty(MinXYZ)==1 , and(length(MinXYZ)==1 , MinXYZ == i)==1 )==1 % if no matches, or matches itself (very short fiber)   ======>>>> UNIQUE VERIFIED  
            % do nothing
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
        elseif length(MinXYZ)==1 % if single match                                                                                          ======>>>> UNIQUE VERIFIED
            % If the max threshold between the in-plane angles and the out-plane angles is 10 degrees 
            if and( abs(DataOut1(i,7)-DataOut1(MinXYZ,7))<fib_stitch_thresh , abs(DataOut1(i,8)-DataOut1(MinXYZ,8))<fib_stitch_thresh  )==1 %check angles match 
                
                % make sure the last slice is always greater than the candidate
                if DataOut1(i,6) < DataOut1(MinXYZ,3)  % z component of the current must always be less than z component of the next                
                StitchMap(j,:)=[fibers(i) fibers(MinXYZ)]; %store the fiber indices that should be stitched
                j=j+1; %count
                end  
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
            % If the fiber index is itself part of the repeating index then exclude the repeating index                                      ======>>>> UNIQUE VERIFIED     
        elseif length(setdiff(MinXYZ, i))==1                               % if 2 match, but 1 is itself (short fiber that needs stitching)
            candidate = setdiff(MinXYZ, i);                                % isolate the candiate
            if and( abs(DataOut1(i,7)-DataOut1(candidate,7))<fib_stitch_thresh , abs(DataOut1(i,8)-DataOut1(candidate,8))<fib_stitch_thresh  )==1 %check angles 
                % make sure the last slice is always greater than the candidate
                
                    if DataOut1(i,6) < DataOut1(candidate,3)  % z component of the current must always be less than z component of the next
                    StitchMap(j,:)=[fibers(i) fibers(candidate)]; %stitch
                    j=j+1; %count
                    end      
            end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
        elseif length(MinXYZ)>1                                                    % bracket sign must exist outside the bracket 
            candidates = setdiff(MinXYZ, i);                                       % isolate the candiate
            compAng = abs(bsxfun(@minus,DataOut1(i,7:8),DataOut1(candidates,7:8)));   %  difference in angles
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%             
            % if the same fiber yields a minimum angle and the differences are less than 10
            if and( and(find(compAng(:,1)==min(compAng(:,1)))==find(compAng(:,2) == min(compAng(:,2))), min(compAng(:,1))<fib_stitch_thresh) ==1 , min(compAng(:,2))<fib_stitch_thresh ) ==1 
                candidate=candidates( find(compAng(:,1)==min(compAng(:,1))) );     % isolate the final candidate

                % make sure the last slice is always greater than the candidate (MAKE SINGLE) 
                    if any(DataOut1(i,6) < DataOut1(candidate,3))                        % z component of the current must always be less than z component of the next
                       candidate = candidate(DataOut1(i,6) < DataOut1(candidate,3));     % extract only those that satisfy the above 
                       [~,min_dist_idx] = min(sqrt(diff(candidate,1).^2 + diff(candidate,2).^2 + diff(candidate,3).^2)); %                   ======>>>> UNIQUE VERIFIED
                        StitchMap(j,:)=[fibers(i) fibers(candidate(min_dist_idx))];      % stitch
                        j=j+1; %count
                    end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                 
            % if the minimum in both cases indicate different fibers and all are within the minimum distance of 10 degrees (Case A)
            elseif and(find(compAng(:,1)== min(compAng(:,1)))~=find(compAng(:,2)==min(compAng(:,2)))  ,  all((reshape(compAng,[],1)<=fib_stitch_thresh)==1))
                
                [~,min_dist_idx] = min(sqrt(diff(candidates,1).^2 + diff(candidates,2).^2 + diff(candidates,3).^2));                      %  ======>>>> UNIQUE VERIFIED
                candidate = candidates(min_dist_idx); 
                
                % make sure the last slice is always greater than the candidate
                    if DataOut1(i,6) < DataOut1(candidate,3)       %  z component of the current must always be less than z component of the next
                    StitchMap(j,:)=[fibers(i) fibers(candidate)];  %  stitch
                    j=j+1; % count
                    end  
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%          
            % if the minimum in both cases indicate different fibers but one of them has values all are within the minimum distance of 10 degrees (Case B)
            elseif and(find(compAng(:,1)== min(compAng(:,1)))~=find(compAng(:,2)==min(compAng(:,2)))  ,  any(all(compAng<=fib_stitch_thresh,2) == 1)) 
             
            % extract candidates whose angles satisfies both the in-plane and out of plane requirement
                candidate = candidates(all(compAng<=fib_stitch_thresh,2)); 
              
                if numel(candidate) > 1     
                [~,idx] = min(min(compAng(all(compAng<=fib_stitch_thresh,2),:),[],2));                                                   %  ======>>>> UNIQUE VERIFIED
                candidate = candidate(idx);
                end
    
                % make sure the last slice is always greater than the candidate
                
                    if DataOut1(i,6) < DataOut1(candidate,3)  % z component of the current must always be less than z component of the next
                    StitchMap(j,:)=[fibers(i) fibers(candidate)]; %stitch
                    j=j+1; %count
                    end               
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              
            % if one candidate yields minimum in both angles but each of the minimum does not satisfy the threshold (Case C)
            % Check the length
            elseif find(compAng(:,1)== min(compAng(:,1)))==find(compAng(:,2)==min(compAng(:,2)))
                candidate=candidates(compAng(:,1)== min(compAng(:,1)));             % extract the  one candidate that yields minimum in both angles but each of the minimum does not satisfy the threshold

%                 % Check to see the length of the fibers 
%                 if norm(DataOut1(candidate,1:3) - DataOut1(candidate,4:6)) <= fib_stitch_thresh_length
%                 % make sure the last slice is always greater than the candidate                
%                     if DataOut1(i,6) < DataOut1(candidate,3)  % z component of the current must always be less than z component of the next
%                     StitchMap(j,:)=[fibers(i) fibers(candidate)]; %stitch
%                     j=j+1; %count
%                     end 
%                 end 
                
                A = DataOut1(candidate,1:3);
                B = DataOut1(candidate,4:6);
                
                if any(sqrt(sum((A-B).^2,2)) <= fib_stitch_thresh_length)
                         candidate = candidate(sqrt(sum((A-B).^2,2)) <= fib_stitch_thresh_length);  
                      if any(DataOut1(i,6) < DataOut1(candidate,3))                        % z component of the current must always be less than z component of the next
                               candidate = candidate(DataOut1(i,6) < DataOut1(candidate,3));     % extract only those that satisfy the above 
                               [~,min_dist_idx] = min(sqrt(diff(candidate,1).^2 + diff(candidate,2).^2 + diff(candidate,3).^2)); %                   ======>>>> UNIQUE VERIFIED
                               StitchMap(j,:)=[fibers(i) fibers(candidate(min_dist_idx))];      % stitch
                               j=j+1; %count
                      end
                end
                     
  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%               

            else                
%                 issue_case(i) = 1;               
% %               This part needs work:
%                 disp('Help! I have multiple matches and none have a good angle match!')
%                 disp(i)
%                 disp(MinXYZ)
%                 disp(compAng)
            end
            
        end
    end
%     disp(['I stitched ', num2str(j-1), ' cases']);
%     toc
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% StitchMap_TMP = StitchMap;
% Extract the common indices and lump them as one fiber 
% (a) First conocatenate 
% (b) Use while loop for the to extract the indices of the fibers 

idx_del          =   zeros(1,size(StitchMap,1));
sub_count        =   1;
cat_vec          =   cell(1,size(StitchMap,1));   
tic 

while ~isempty(StitchMap)     
%     sprintf('%s%d','remaining rows of stitch map is ',size(StitchMap,1))
    
   idx_del_tmp   =   idx_del;
   loop_f        =   1;
   idx_end       =   StitchMap(1,2);   
   table_var     =   zeros(1,2);             %  initialize the table_var 
    
   table_var(1,:)    =   StitchMap(1,:);
   idx_del_tmp(1)    =   1;
   count_2           =   1;
    
    while loop_f == 1  
        count_2 = count_2 + 1;  
        %   Trailing end of fiber does not have to be stitched 
        if isempty(StitchMap(StitchMap(:,1) == idx_end,:))     
            loop_f = 0; 
        else             
            %   Trailing end of fiber has to be stitched
            idx_del_tmp(count_2) = find(StitchMap(:,1) == idx_end);
            table_var(count_2,:) = StitchMap(StitchMap(:,1) == idx_end,:);
            %    Update the idx_end variable
            idx_end = table_var(end); 
            
            % prevent infinite loop here 
            tab_reshape = reshape(table_var,1,[]);
            if ismember(idx_end, tab_reshape(1:end-1))
                loop_f = 0;
            end
            
        end                  
    end 
    
    idx_del_tmp(idx_del_tmp == 0)     =     [];
    cat_vec{sub_count}                =     unique(reshape(table_var,1,[]));
    StitchMap(idx_del_tmp,:)          =     [];
    sub_count                         =     sub_count + 1;    
end
toc 

% Remove the empty cells 
cat_vec(cell2mat(cellfun(@(x)isempty(x),cat_vec,'uni',0)))  =   [];

%% OUTPUT 1: IDX_to_keep
disp('Begin IDX_to_keep phase')
% Create a variable that takes into account the indices of the fibers after stitching 
IDX_to_keep = [(num2cell(fibers(~ismember(fibers,[cat_vec{:}])),2))'     cat_vec]; 
disp('End IDX_to_keep phase')

%% OUTPUT 2: min_slice_Amalgam
disp('Begin min_slice_Amalgam phase')
z                     =     cellfun(@(x)[x(1) x(end)],ze,'uni',0);
min_slice_AMALGAM     =     zeros(1,size(IDX_to_keep,2));

for ii = 1:size(IDX_to_keep,2)    
    if  numel(IDX_to_keep{ii})>1
        tmp = z(IDX_to_keep{ii});
        tmp = [tmp{:}];
        min_slice_AMALGAM(ii)    =     min(tmp);
    else 
        min_slice_AMALGAM(ii)    =     min(z{IDX_to_keep{ii}});
    end
end
disp('End min_slice_Amalgam phase')

%% OUTPUT 3: FINAL_CENTROID_AMALGAM, fib_diam_AMALGAM  and  Linear_Index_AMALGAM
% For the FINAL_CENTROID_AMALGAM interpolate through the fibers as a stitch
Centroid_MID_cell_LUMP_ORIG_UP         =     cellfun(@(x,y)[x(:,1) x(:,2) (y(1):y(end))'],Centroid_MID_cell_LUMP_ORIG,z,'uni',0);
FINAL_CENTROID_AMALGAM                 =     cell(1,size(IDX_to_keep,2));
Linear_Index_AMALGAM                   =     cell(1,size(IDX_to_keep,2));
fib_diam_AMALGAM                       =     zeros(1,size(IDX_to_keep,2));


disp('Begin FINAL_CENTROID_AMALGAM, fib_diam_AMALGAM  &  Linear_Index_AMALGAM phase')
for ii = 1:size(IDX_to_keep,2) 
    
    if numel(IDX_to_keep{ii})>1                        
        Cent_VEC       =   Centroid_MID_cell_LUMP_ORIG_UP(IDX_to_keep{ii});
        Cent_VEC       =   vertcat(Cent_VEC{:});    
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fin_diam_stitched
        diam_vec                 =    minRadius_idxval_LUMP_cat(IDX_to_keep{ii});
        diam_vec                 =    vertcat(diam_vec{:}); 
        std_diam                 =    std(diam_vec);
        mean_diam                =    mean(diam_vec);        
        fib_diam_AMALGAM(ii)     =    mean(diam_vec(diam_vec>=(mean_diam-std_diam) & diam_vec<=(std_diam+mean_diam)));
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Linear_Index_AMALGAM
        Linear_Index_AMALGAM_tmp         =    Linear_Index_center_GLOBAL_LUMP(IDX_to_keep{ii});
        Linear_Index_AMALGAM{ii}         =    vertcat(Linear_Index_AMALGAM_tmp{:});
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  
        [~,idx]                          =       sort(Cent_VEC(:,3));       % sort just the third column
        Cent_VEC                         =       Cent_VEC(idx,:);           % sort the whole matrix using the sort indices    
        Cent_VEC_tmp                     =       Cent_VEC;
        
        % Ensure that the query_vec_points always start from the least to the max 
        query_vec_points                 =       min(Cent_VEC_tmp(:,3)) : max(Cent_VEC_tmp(:,3));
    
        clearvars idx_val 
        [idx_val(:,2),idx_val(:,1)]      =       hist(Cent_VEC_tmp(:,3),unique(Cent_VEC_tmp(:,3)));     
        Cent_VEC_tmp                     =       Cent_VEC_tmp(ismember(Cent_VEC_tmp(:,3),idx_val(idx_val(:,2)==1,1)),:);    
     
        % Ensures that the interpolation of the centroids is always fixed.
            if size(Cent_VEC_tmp,1) == 1 
                Cent_VEC(2:end-1,:) = [];
            else    
                Cent_VEC = Cent_VEC_tmp;   
            end
        
        sample_points_cr   =  Cent_VEC(:,3);
        sample_values_cr   =  Cent_VEC(:,1:2);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Perform interpolation between the centroids
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        clearvars result
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % interpolation for centroids columns
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        result(:,1)      =    interp1(sample_points_cr,sample_values_cr(:,1),query_vec_points','linear','extrap');
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % interpolation for centroids rows
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        result(:,2)                      =    interp1(sample_points_cr,sample_values_cr(:,2),query_vec_points','linear','extrap');
        FINAL_CENTROID_AMALGAM{ii}       =    [result   query_vec_points'];
                     
    else
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % FINAL_CENTROID_AMALGAM
        FINAL_CENTROID_AMALGAM{ii}       =    Centroid_MID_cell_LUMP_ORIG_UP{IDX_to_keep{ii}};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Linear_Index_AMALGAM
        Linear_Index_AMALGAM{ii}         =    Linear_Index_center_GLOBAL_LUMP{IDX_to_keep{ii}};
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % fin_diam_stitched
        diam_vec                         =    minRadius_idxval_LUMP_cat{IDX_to_keep{ii}};
        std_diam                         =    std(diam_vec);
        mean_diam                        =    mean(diam_vec);
        fib_diam_AMALGAM(ii)             =    mean(diam_vec(diam_vec>=(mean_diam-std_diam) & diam_vec<=(std_diam+mean_diam)));

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%        
    end
             
end

disp('End FINAL_CENTROID_AMALGAM, fib_diam_AMALGAM  &  Linear_Index_AMALGAM phase')

varargout{1}   =   IDX_to_keep;
varargout{2}   =   fib_diam_AMALGAM;
varargout{3}   =   Linear_Index_AMALGAM;
varargout{4}   =   FINAL_CENTROID_AMALGAM;
varargout{5}   =   min_slice_AMALGAM;

end


%%
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%            
%                  
%             [xx,yy,zz] = cellfun(@(x)ind2sub(size_length,Linear_Index_center_GLOBAL_LUMP{x}),num2cell([candidate;fibers(i)],2),'uni',0);
%             
%             max(vertcat(xx{:}))
%             max(vertcat(yy{:}))   
%             max(vertcat(zz{:}))  
%             
%             min(vertcat(xx{:}))
%             min(vertcat(yy{:}))   
%             min(vertcat(zz{:})) 
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  



%% OUTPUT 4: Linear_Index_Amalgam
% % Linear_Index_Amalgam also easy to do 
% Linear_Index_AMALGAM  =    Linear_Index_center_GLOBAL_LUMP(IDX_to_keep);
% Linear_Index_AMALGAM  =    vertcat(Linear_Index_AMALGAM{:});

%% OUTPUT 5: 


%%
% ========================= % 
% Terms for the variables % 
% ========================= % 

% IDX_to_keep;
% fib_diam_AMALGAM;
% Linear_Index_AMALGAM;
% FINAL_CENTROID_AMALGAM;
% min_slice_AMALGAM;
% duplicate_vector;


%%

%         % Extract the slice heights of elements within each cell 
%         slice_height_info = cellfun(@(x)x(:,3),Centroid_MID_cell_LUMP_ORIG_UP(IDX_to_keep{ii}),'uni',0);
%         
%         % Run a loop that zeros indices of common slice heights
%         for kk = 1:size(slice_height_info,2)-1
%             common_members = slice_height_info{kk}(ismember(slice_height_info{kk},slice_height_info{kk+1}));
%             slice_height_info{kk}(ismember(slice_height_info{kk},common_members)) = 0;
%             slice_height_info{kk+1}(ismember(slice_height_info{kk+1},common_members)) = 0;
%         end
%         
%         % Extract the remaining slice indices after the zeroing and use that as interpolants to perform the linear interpolation 
%         slice_height_info_tmp = vertcat(slice_height_info{:});
%         slice_height_info_tmp(slice_height_info_tmp==0) = [];



%%
% parfor ult_count = 1:numel(avail_terms)
%     
%    sprintf('%s%d','ult_count_',ult_count)
%    
% %% STEP (II) FORWARD SEARCH ANALYSIS:    
% % clearvars   FINAL_NEW_CENTROID_F         keep_idx_f         negatives_extract_f
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % for degbugging purposes
% % if ult_count == 1188 % This is a scenario of two fibers are to be connected to one large stock fiber 
% %     break
% % end
%  
% % if ult_count == 1787
% %     break
% % end
% 
% % if avail_terms(1) == 752
% %     break
% % end
% 
% % if min(z_TMP{895}) ~= 73
% %     break
% % end
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% %% STEP (III) Set TEMPORARY variables that keep original dimensions
% % Set the z variable to extract the extremities of the slice locations 
% z                                      =     cellfun(@(x)[x(1) x(end)],ze,'uni',0);
% z_TMP                                  =     z;
% 
% % Update the centroid of the constituent ellipses using the appropriate heights 
% Centroid_MID_cell_LUMP_ORIG_UP         =     cellfun(@(x,y)[x(:,1) x(:,2) (y(1):y(end))'],Centroid_MID_cell_LUMP_ORIG,z,'uni',0);
% Centroid_MID_cell_LUMP_TMP             =     Centroid_MID_cell_LUMP_ORIG_UP;
% STATS_Pix_LOC_LUMP_UP                  =     STATS_Pix_LOC_LUMP;
% minRadius_idxval_LUMP_cat_TMP          =     minRadius_idxval_LUMP_cat;
% Linear_Index_center_GLOBAL_LUMP_TMP    =     Linear_Index_center_GLOBAL_LUMP;
% 
% % DESCRIPTION OF THE FORWARD STITCH ANALYSIS:
% % This framework is such that for each fiber;
% % (a)  Project through the region above to find if a fiber exists or if a rogue region exists or if an empty space exists. 
%        % (a - i)    If a FIBER EXIST: makes sure that invoke the sieve mechanism on the basis of original connected pairs deletion, in_plane distance thresh and orientation thresh. 
%        % (a - ii)   If a ROGUE REGION EXISTS: progressively move through the volume till you find a fiber (Take note of the ends of the fiber-rogue_region-fiber combination).
%        % (a - iii)  If an EMPTY REGION EXISTS: project through 3 slices in the dxn of the line fit to see of a fiber or rogue region exists.      
% % (b)  For the in-plane analysis perform a bi-driectional analysis to resolve the dichotomy instances.
% % (e)  Extract the corresponding centroid positions into an appropriate variable. 
% 
% % if isempty(horzcat(IDX_to_keep{:}),avail_terms(ult_count))
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% c_f                             =    avail_terms(ult_count);               % Make change here to loop through
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % BEGIN "IF-STATEMENT" HERE:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% % if isempty(ismember(IDX_to_keep_TMP,avail_terms(ult_count)))
% 
% radius                          =    4;
% MAIN_counter                    =    1;
% keep_idx_f                      =    cell(1,10);                           % Chnge Here
% FINAL_NEW_CENTROID_F            =    cell(1,10);
% negatives_extract_f             =    cell(1,10);
% keep_idx_f{MAIN_counter}        =    c_f;
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        
%                       % Check to see if the first element of the variable avail_term is non-empty 
%                       % NOTE: empty means no fiber stitches to the current fiber 
%                       
%                        while ~isempty(c_f)      
%                        
%                             if z_TMP{c_f}(end) < roof_limit    % Check if the end of the seed fiber flushes with the ceiling of the 3D volume 
%                                 
%                                 % The function CENTERS_EXTRAPOLATE_LOOP_FUNC_NEW_UPDATE that scans through the volume of the TEMP_MAT_NEW and extracts the index of a fiber in the volume 
%                                 % (a) Function stops upon the first scan when it extrats fiber/s 
%                                 % (b) Function projects through 3slices for initial empty scan 
%                                 % (c) Function projects through the rogue volume to extract either fiber/s or it ends within the rogue volume  
%                                 
%                                  zf = z_TMP{c_f};
%                                  halt_var = 1;
%                                  slice_limit = 0;
%                                  query_path = 2;
%                                  
%                                  [circlePixels_points_IDX,main_res_extract,main_res,...
%                                   query_var,query_var_tmp,crude_bin_mat]              =        CENTERS_EXTRAPOLATE_LOOP_FUNC_NEW_UPDATE(halt_var,zf,slice_limit,c_f,radius,columnsInImage,...
%                                                                                                                  rowsInImage,size_length,Centroid_MID_cell_LUMP_TMP,query_path,TEMP_MAT_NEW,Structure_Rogue,z_TMP,ORIGINAL_CRUDE_BIN_MAT); 
%                                     
%                                     % ==> Check if potential fibers are there to be stitched ?                                                                         
%                                     if isempty(query_var)
%                                             
%                                             % if the potential fibers do not exist then two scenarions may exist 
%                                             % (a) probe fiber may impinge a Rogue region 
%                                             % (b) probe fiber may not impinge a Rogue region
%                                               
%                                             if any(vertcat(query_var_tmp{:})<0)    % If probe fiber impinges a Rogue region 
%                                             
%                                                 % Perform the following:
%                                                 
%                                                 %(a)  Extract the indices of the rogue regions 
%                                                 tmp_cell = cellfun(@(x)x(x<0),query_var_tmp,'uni',0);
%                                                 tmp_cell = abs(vertcat(tmp_cell{:})');
% 
%                                                 %(b)  Extract the maximum slice height of the rogue region(s)                                            
%                                                 [max_rogue,~] = max(cell2mat(cellfun(@(x)max(x),Structure_Rogue.rogue_ends(tmp_cell),'uni',0)));                                            
%                                                  main_res(main_res(:,3) > max_rogue,:) = [];
% 
%                                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                           
%                                                 %%%%%% PROJECTING THROUGH THE ROGUE VOLUME IF NEED BE %%%%%%%
% 
%                                                 % For original centroids with size less than the max of the projected region perform an extrapolation
%                                                         % Note here that we invoke MINOR_EXTRAPOLATION_FUNC if original centroid has more than one point else use updated centroids 
%                                                         if size(main_res,1)>size(Centroid_MID_cell_LUMP_ORIG_UP{c_f},1)
%                                                                if size(Centroid_MID_cell_LUMP_ORIG_UP{c_f},1)>1
%                                                             NEW_CENTROID =  MINOR_EXTRAPOLATION_FUNC(main_res,Centroid_MID_cell_LUMP_ORIG_UP{c_f});
%                                                                elseif  size(Centroid_MID_cell_LUMP_ORIG_UP{c_f},1)==1
%                                                             NEW_CENTROID =   main_res;     
%                                                                end
%                                                 % For original centroids with size equal to the max of the projected region do not perform an extrapolation 
%                                                         else 
%                                                             NEW_CENTROID   = Centroid_MID_cell_LUMP_ORIG_UP{c_f};   
%                                                         end
% 
%                                                         c_f             =    query_var;
%                                                 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%                                             
%                                             % Indicates Rogue regions do not exist                                            
%                                             elseif ~any(vertcat(query_var_tmp{:})<0)
%                                             % This suggests that no rogue region was found for this condition to be true     
%                                             % (NOTE: main_res is the  projected centroid a slice beyond the fiber or beyond a rogue region attached to the fiber) 
%                                                                    
%                                                     NEW_CENTROID    =     Centroid_MID_cell_LUMP_ORIG_UP{c_f};
%                                                     c_f             =     query_var;
%                                             end
%                                             
%                                     elseif ~isempty(query_var)
%                                       
%                                             %%%%%%% CHECK 1: SLICE LEVEL CHECK %%%%%%%%
%                                             % Delete all fibers whose root ellipses begin on slices lower than the projected slice 
%                                             % Note these are fibers that were either connected pairs or fibers that touch at some junction
%                                             % the fibers to be stitched are the ones that always begin at the projected slice ar some slices higher due to the rogue regions
%                                            
%                                            %% Verify that the fibers 
%                                            % (a) Do not encase other fibers and are not due to oversegmentation of an in-plane ellipses
%                                            % (b) Are not a result of oversegmented 3D fibers 
%                                            
%                                          [query_var,query_var_over_segTMP,d_vec_angles_overTMP]   =   ST_query_over_segTMP_f(c_f,query_var,Crude_Idx_Split_Idx_var,z_TMP,ze,...
%                                                                                                                           Linear_Index_center_GLOBAL_LUMP_TMP,STATS_Pix_LOC_LUMP_UP,...
%                                                                                                                           direction_vectors,size_length);
%                                           
%                                           %% Here we shall look into the limit of the fibers' ends 
%                                            % Note that mes_res_extract is projected centroids by some distance.
%                                            % There is therefore the need to correct for the extraoneous projection of the crude bin fiber track. 
%                                             
%                                           [query_var,~,main_res]                   =     ST_query_further_analysis_f(query_var,crude_bin_mat,Linear_Index_center_GLOBAL_LUMP_TMP,size_length,...
%                                                                                                                                     main_res,main_res_extract,rowsInImage,columnsInImage,radius,z_TMP,...
%                                                                                                                                     Centroid_MID_cell_LUMP_TMP,TEMP_MAT_NEW,Structure_Rogue,ORIGINAL_CRUDE_BIN_MAT,...
%                                                                                                                                     direction_vectors,c_f,inplane_thresh,ang_thresh,circlePixels_points_IDX);
%                                        
% %% Compare the fibers for the shared volume and the unshared volume 
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % 
%                                         [NEW_CENTROID,c_f,z,z_TMP,Centroid_MID_cell_LUMP_ORIG_UP,...
%                                          Centroid_MID_cell_LUMP_TMP,STATS_Pix_LOC_LUMP_UP]     =    ST_choose_over_or_query(query_var_over_segTMP,query_var,direction_vectors,d_vec_angles_overTMP,...
%                                                                                                                                                      Centroid_MID_cell_LUMP_ORIG_UP,Centroid_MID_cell_LUMP_TMP,z_TMP,z,c_f,...
%                                                                                                                                                      STATS_Pix_LOC_LUMP_UP,ang_thresh,main_res);
% 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                         
%                                     end
%                                     
%                                 % ADD TO RESULT TO THE LIST FOR THAT PROJECTION                                  
%                                 FINAL_NEW_CENTROID_F{MAIN_counter}           =     NEW_CENTROID;
%                                 tmp_cell                                     =     cellfun(@(x)x(x<0),query_var_tmp,'uni',0);
%                                 negatives_extract_f{MAIN_counter}            =     [tmp_cell{:}];
%                                 
%                             else   % Check if the end does not flush with the ceiling of the 3D volume 
%                                 
%                                 negatives_extract_f{MAIN_counter}            =     [];
%                                 FINAL_NEW_CENTROID_F{MAIN_counter}           =     Centroid_MID_cell_LUMP_ORIG_UP{c_f};
%                                 keep_idx_f{MAIN_counter}                     =     c_f;
%                                 c_f                                          =     [];
%                             end
%                             
%                                 MAIN_counter                                 =     MAIN_counter + 1;                       
%                                 keep_idx_f{MAIN_counter}                     =     c_f;
% %                               clearvars probe_var_tmp     probe_var     tmp_cell
%                        end         
% 
%            %% BACKWARD SEARCH           
%             c_f = avail_terms(ult_count);
% %           clearvars FINAL_NEW_CENTROID_B    keep_idx_b    negatives_extract_b
% 
%             MAIN_counter                           =    1;
%             
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             keep_idx_b                      =    cell(1,10);                %   Chnge Here  
%             FINAL_NEW_CENTROID_B            =    cell(1,10);
%             negatives_extract_b             =    cell(1,10);
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                   
%                       
%             keep_idx_b{MAIN_counter}               =    c_f;
%             
%                        % Check to see if the first element of the variable avail_term is non-empty 
%                        while ~isempty(c_f)
% 
%                             if z_TMP{c_f}(1)  >  1   % Check if the beginning slice of the seed fiber flushes with the bottom of the 3D volume 
%                                 
%                                 % Create a function that scans through the volume iteratively till it meets another fiber or sets of fibers and then it stops.
%                                 % In other words the script scans through rogue volumes till it meets another fiber and then it stops
%                                  zf            =   z_TMP{c_f};
%                                  halt_var      =   1;
%                                  slice_limit   =   0; 
%                                  query_path    =   1; 
%                                                                
%                                  [circlePixels_points_IDX,main_res_extract,main_res,...
%                                   query_var,query_var_tmp,crude_bin_mat]       =       CENTERS_EXTRAPOLATE_LOOP_FUNC_NEW_UPDATE(halt_var,zf,slice_limit,c_f,radius,columnsInImage,...
%                                                                                                  rowsInImage,size_length,Centroid_MID_cell_LUMP_TMP,query_path,TEMP_MAT_NEW,Structure_Rogue,z_TMP,ORIGINAL_CRUDE_BIN_MAT);    
%                                                                                                       
%                                     if isempty(query_var)
%                                     
%                                             % Indicates Rogue regions exist
%                                             if any(vertcat(query_var_tmp{:})<0)                                           
%                                             
%                                             % Extract the indices of the rogue regions
%                                             tmp_cell = cellfun(@(x)x(x<0),query_var_tmp,'uni',0);
%                                             tmp_cell = abs(vertcat(tmp_cell{:})');
%                                             
%                                             % Extract the maximum slice height of the rogue region(s)
%                                             [min_rogue,~] = min(cell2mat(cellfun(@(x)min(x),Structure_Rogue.rogue_ends(tmp_cell),'uni',0))); 
%                                             main_res(main_res(:,3) < min_rogue,:) = [];
%                                             
%                                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                           
%                                             %%%%%% PROJECTING THROUGH THE ROGUE REGION IF NEED BE %%%%%%%
%                                             % For original centroids with size less than the max of the projected region perform an extrapolation
%                                             
%                                                     % Note here that we invoke MINOR_EXTRAPOLATION_FUNC if original centroid has ore than one point else use updated centroids 
%                                                     if size(main_res,1)>size(Centroid_MID_cell_LUMP_ORIG_UP{c_f},1)
%                                                         if size(Centroid_MID_cell_LUMP_ORIG_UP{c_f},1) > 1
%                                                     NEW_CENTROID =  MINOR_EXTRAPOLATION_FUNC(main_res,Centroid_MID_cell_LUMP_ORIG_UP{c_f});
%                                                         elseif size(Centroid_MID_cell_LUMP_ORIG_UP{c_f},1) == 1
%                                                     NEW_CENTROID =  main_res;
%                                                         end
%                                                         
%                                                     % For original centroids with size equal to the max of the projected region do not perform an extrapolation     
%                                                     else 
%                                                     NEW_CENTROID    = Centroid_MID_cell_LUMP_ORIG_UP{c_f};   
%                                                     end
%                                                     
%                                                     c_f             =    query_var; 
%                                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                             
%                                             
%                                             % Indicates Rogue regions do not exist 
%                                             elseif ~any(vertcat(query_var_tmp{:})<0)
%                                             % This suggests that no rogue region was found for this condition to be true     
%                                             % (NOTE: main_res is the  projected centroid a slice beyond the fiber or beyond a rogue region attached to the fiber) 
% 
%                                                     NEW_CENTROID    =     Centroid_MID_cell_LUMP_ORIG_UP{c_f};
%                                                     c_f             =     query_var;                                            
% %                                             % An empty probe_var indicates that no region of fiber was found after the projection
% %                                             % (NOTE: main_res is the  projected centroid a slice beyond the fiber or beyond a rogue region attached to the fiber)                                            
% %                                             main_res(main_res(:,3) < z_TMP{c_f}(1),:) = [];
% % %                                             NEW_CENTROID = main_res;
% %                                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                           
% %                                             %%%%%% KEEPING THE ORIGINAL CENTROID %%%%%%%
% %                                                     if size(main_res,1)>size(Centroid_MID_cell_LUMP_ORIG_UP{c_f},1)
% %                                                     if size(Centroid_MID_cell_LUMP_ORIG_UP{c_f},1) > 1
% %                                                     NEW_CENTROID =  MINOR_EXTRAPOLATION_FUNC(main_res,Centroid_MID_cell_LUMP_ORIG_UP{c_f});
% %                                                     elseif size(Centroid_MID_cell_LUMP_ORIG_UP{c_f},1) == 1
% %                                                     NEW_CENTROID = main_res;
% %                                                     end
% %                                                     else 
% %                                                     NEW_CENTROID =  Centroid_MID_cell_LUMP_ORIG_UP{c_f};   
% % 
% %                                                     end
%                                             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                                                                                           
%                                             end
%                                             
%                                     elseif ~isempty(query_var) 
%                                             %%%%%%% CHECK 1: SLICE LEVEL CHECK %%%%%%%%
%                                             % Delete all fibers whose root ellipses begin on slices lower than the projected slice 
%                                             % Note these are fibers that were either connected pairs or fibers that touch at some junction
%                                             % the fibers to be stitched are the ones that always begin at the projected slice ar some slices higher due to the rogue regions
% 
%                                             % Verify that the fibers 
%                                             % (a) Do not encase other fibers and the not due oversegmentation of an in-plane ellipses
%                                             
%                                          [query_var,query_var_over_segTMP,d_vec_angles_overTMP]   =   ST_query_over_segTMP_b(c_f,query_var,Crude_Idx_Split_Idx_var,z_TMP,ze,...
%                                                                                                                           Linear_Index_center_GLOBAL_LUMP_TMP,STATS_Pix_LOC_LUMP_UP,...
%                                                                                                                           direction_vectors,size_length);
%                                                                                
%                                           %% Here we shall look into the limit of the fibers' ends 
%                                            % Note that mes_res_extract is projected centroids by some distance.
%                                            % There is therefore the need to correct for the extraoneous projection of the crude bin fiber track. 
% 
%                                          [query_var,~,main_res]                   =     ST_query_further_analysis_b(query_var,crude_bin_mat,Linear_Index_center_GLOBAL_LUMP_TMP,size_length,...
%                                                                                                                                    main_res,main_res_extract,rowsInImage,columnsInImage,radius,z_TMP,...
%                                                                                                                                    Centroid_MID_cell_LUMP_TMP,TEMP_MAT_NEW,Structure_Rogue,ORIGINAL_CRUDE_BIN_MAT,...
%                                                                                                                                    direction_vectors,c_f,inplane_thresh,ang_thresh,circlePixels_points_IDX);
%                                             
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
% %% Compare the fibers for the shared volume and the unshared volume 
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % if both are not empty select the best 
% % if probe is not empty and O is empty 
% % if probe is empty and O is not empty 
% % if probe and O are not empty
% 
%                                     [NEW_CENTROID,c_f,z,z_TMP,Centroid_MID_cell_LUMP_ORIG_UP,...
%                                      Centroid_MID_cell_LUMP_TMP,STATS_Pix_LOC_LUMP_UP]     =    ST_choose_over_or_query(query_var_over_segTMP,query_var,direction_vectors,d_vec_angles_overTMP,...
%                                                                                                                         Centroid_MID_cell_LUMP_ORIG_UP,Centroid_MID_cell_LUMP_TMP,z_TMP,z,c_f,...
%                                                                                                                         STATS_Pix_LOC_LUMP_UP,ang_thresh,main_res);
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                     end
%                                     
%                                 % ADD TO RESULT TO THE LIST FOR THAT PROJECTION                                     
%                                 FINAL_NEW_CENTROID_B{MAIN_counter} = NEW_CENTROID;
%                                 tmp_cell = cellfun(@(x)x(x<0),query_var_tmp,'uni',0);
%                                 negatives_extract_b{MAIN_counter}    = [tmp_cell{:}];
%                                 
%                             else  % Check if the beginning of the slice fiber does not flush with the ceiling of the 3D volume
% 
%                                 negatives_extract_b{MAIN_counter}    = [];
%                                 FINAL_NEW_CENTROID_B{MAIN_counter}           =     Centroid_MID_cell_LUMP_ORIG_UP{c_f};
%                                 keep_idx_b{MAIN_counter} = c_f;
%                                 c_f = [];
%                             end
% 
%                                 MAIN_counter = MAIN_counter + 1;
%                                 keep_idx_b{MAIN_counter} = c_f;                                
% %                               clearvars     probe_var_tmp      probe_var      tmp_cell     
%                        end       
% 
%             %% LUMP TOGETHER 
%             FINAL_NEW_CENTROID_B                   =     vertcat(FINAL_NEW_CENTROID_B{:});
%             [~,idxb]                               =     sort(FINAL_NEW_CENTROID_B(:,end),'ascend');
%             FINAL_NEW_CENTROID_B                   =     FINAL_NEW_CENTROID_B(idxb,:);  
%             FINAL_NEW_CENTROID_F                   =     vertcat(FINAL_NEW_CENTROID_F{:});
%             [~,idxf]                               =     sort(FINAL_NEW_CENTROID_F(:,end),'ascend');
%             FINAL_NEW_CENTROID_F                   =     FINAL_NEW_CENTROID_F(idxf,:);
%             
%             % Delete the fibers that repeat in the variable FINAL_NEW_CENTROID_B and FINAL_NEW_CENTROID_F
%             FINAL_NEW_CENTROID_F(ismember(FINAL_NEW_CENTROID_F(:,3),FINAL_NEW_CENTROID_B(:,3)),:)   = [];            
%             MAIN_FINAL_NEW_CENTROID                =      vertcat(FINAL_NEW_CENTROID_B,FINAL_NEW_CENTROID_F);                                      
%             MAIN_FINAL_NEW_CENTROID(:,3)           =      (1:size(MAIN_FINAL_NEW_CENTROID,1))';
%             
%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               % % Ensure all data is in row form
%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              [keep_idx_f]             =     Horz_func(keep_idx_f);
%              [keep_idx_b]             =     Horz_func(keep_idx_b);
%              [negatives_extract_f]    =     Horz_func(negatives_extract_f);
%              [negatives_extract_b]    =     Horz_func(negatives_extract_b);            
%              keep_final               =     unique([horzcat(keep_idx_f{:})              horzcat(keep_idx_b{:})]);           
%              keep_final_negative      =     unique([horzcat(negatives_extract_b{:})     horzcat(negatives_extract_f{:})]);  
%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%    
%             % Note that if keep_final_negative is empty then it means no rogue regions exist in that analysis            
%             if isempty(keep_final_negative)
%             min_slice    =   min(cell2mat(cellfun(@(x)min(x),z_TMP(keep_final),'uni',0)));
%             else            
%             min_slice    =   min(min(cell2mat(cellfun(@(x)min(x),z_TMP(keep_final),'uni',0))),...
%                                  min(cell2mat(cellfun(@(x)min(x),Structure_Rogue.rogue_ends(abs(keep_final_negative)),'uni',0))));
%             end          
%             IDX_to_keep{ult_count}      =    keep_final;
%             
%             % extract the pixel locations and the minor axis lengths of the sticthed fibers 
%             FINAL_CENTROID_stitched{ult_count}   =   MAIN_FINAL_NEW_CENTROID;           
% %           fib_diam_stitched(ult_count)         =   mean(cell2mat(cellfun(@(x)x(1),minRadius_idxval_LUMP_cat_TMP(keep_final),'uni',0)));    % Wrong Analysis of the
%             fib_diam_stitched(ult_count)         =   mean(cell2mat(cellfun(@(x)mean(x(x>(mean(x)-std(x)) & x<(std(x)+mean(x)))),minRadius_idxval_LUMP_cat_TMP(keep_final),'uni',0)));    % Correct Analysis of the
%             Linear_Index_stitched{ult_count}     =   vertcat(Linear_Index_center_GLOBAL_LUMP_TMP{keep_final});           
%             min_slice_stitched(ult_count)        =   min_slice;             
% 
% %             mean(cell2mat(cellfun(@(x)x(1),minRadius_idxval_LUMP_cat_TMP(781),'uni',0)))
% %             cell2mat(cellfun(@(x)mean(x>(mean(x)-std(x)) & x<(std(x)+mean(x))),minRadius_idxval_LUMP_cat_TMP(keep_final),'uni',0))            
% %             cellfun(@(x)std(x),minRadius_idxval_LUMP_cat_TMP(keep_final),'uni',0)
% %             cellfun(@(x)x(x>(mean(x)-std(x))),minRadius_idxval_LUMP_cat_TMP(keep_final),'uni',0)
% 
% % else 
% %             IDX_to_keep{ult_count}               =    {};
% %             % extract the pixel locations and the minor axis lengths of the sticthed fibers 
% %             FINAL_CENTROID_stitched{ult_count}   =    {};           
% %             fib_diam_stitched(ult_count)         =    {};           
% %             Linear_Index_stitched{ult_count}     =    {};           
% %             min_slice_stitched(ult_count)        =    {};
% %  
% %             
% % end 
% %
% if numel(keep_final)>1
% duplicate_vector(ult_count) = 1;
% else 
% duplicate_vector(ult_count) = 0;
% end
% 
% % IDX_to_keep_TMP = [IDX_to_keep{:}];
% % IDX_to_keep_TMP = [IDX_to_keep{:}];
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % END "IF-STATEMENT" HERE:
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
%             
% end
% 
% % find(cell2mat(cellfun(@(x)any(x == 1246),IDX_to_keep,'uni',0)))
% % find(cell2mat(cellfun(@(x)any(x == 1246),IDX_to_keep_FINAL,'uni',0)))
% 
% %% Remove all the duplicate elements 
% [IDX_to_keep,fib_diam_stitched,min_slice_stitched,...
%  Linear_Index_stitched,FINAL_CENTROID_stitched]                 =        DUPLICATE_CORRECTION_FUNC(IDX_to_keep,fib_diam_stitched,min_slice_stitched,...
%                                                                                                    Linear_Index_stitched,FINAL_CENTROID_stitched);
% %  Include the metadata relating to the bad fiber terms
% z_TMP                    =     cellfun(@(x)[x(1) x(end)],ze,'uni',0);
% 
% %  Delete the zeros in the fib_diam_AMALGAM and the min_slice_stitched
% fib_diam_stitched(fib_diam_stitched == 0)     =   [];
% min_slice_stitched(min_slice_stitched == 0)   =   [];
% 
% if Exclude_bad == 1
%     
% Linear_Index_AMALGAM      =     [Linear_Index_stitched         Linear_Index_center_GLOBAL_LUMP_UP(bad_fibers)];
% % temp_2                    =     cell2mat(cellfun(@(x)x(1),minRadius_idxval_LUMP_cat(bad_fibers),'uni',0));
% temp_2                    =     cell2mat(cellfun(@(x)mean(x(x>(mean(x)-std(x)) & x<(std(x)+mean(x)))),minRadius_idxval_LUMP_cat(bad_fibers),'uni',0));
% fib_diam_AMALGAM          =     [fib_diam_stitched             temp_2];
% FINAL_CENTROID_AMALGAM    =     [FINAL_CENTROID_stitched       Centroid_MID_cell_LUMP_ORIG(bad_fibers)];
% z_start                   =     cell2mat(cellfun(@(x)min(x),z_TMP(bad_fibers),'uni',0));
% min_slice_AMALGAM         =     [min_slice_stitched            z_start];
% 
% else 
% 
% Linear_Index_AMALGAM      =     Linear_Index_stitched;
% fib_diam_AMALGAM          =     fib_diam_stitched;
% FINAL_CENTROID_AMALGAM    =     FINAL_CENTROID_stitched;
% min_slice_AMALGAM         =     min_slice_stitched;
% 
% end
% 
% 
% varargout{1} = IDX_to_keep;
% varargout{2} = fib_diam_AMALGAM;
% varargout{3} = Linear_Index_AMALGAM;
% varargout{4} = FINAL_CENTROID_AMALGAM;
% varargout{5} = min_slice_AMALGAM;
% varargout{6} = duplicate_vector;
% end
