
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Version created on: 9/2/2018
% To avoid the issue of out of memory error we take into account the cell
% version of the image volumes
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [varargout] = EPAD_FUNCTION_CODE_MASTER_VERSION_MUT_NEW(varargin)

curr_vol    =   varargin{1};
succ_vol    =   varargin{2};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This function ensures that the initial fiber indexing process is fast and very reliable 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Use the EPAD Algorithm to ensure that find the potential fibers 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The implmentation christened "EPAD", we intend to use the arithmetic operations of 
% (a) Exponentiation
% (b) Premultiplication 
% (c) Addition 
% (d) Division 

% Find the overlap of the fiber indices
disp('Find the overlap of the fiber indices') 
exponent = 10^(numel(num2str(max(succ_vol(:)))) + 2);

% Find the sum of the product 
% curr_vol_pre       =   single(curr_vol*exponent);
disp('Find the sum of the product')
curr_succ_SUM      =   (curr_vol*exponent) + succ_vol;

%% This helps bring out the overlap
% Use MATLAB's hist to count the numebr of overlapping indices in the volume 
tic
% Find all the nonzero numbers and their respective occurences 
clearvars p_idx_mod  p_idx_mod_1  curr_succ_overlap
disp('Find all the nonzero numbers and their respective occurences')
[p_idx_mod(:,2),p_idx_mod(:,1)] = hist(nonzeros(curr_succ_SUM(:)),unique(nonzeros(curr_succ_SUM(:))));
toc

%%%%%%%%%%%%%%%%%%%%%%%%
disp('relieve memory == >> curr_succ_SUM')
clearvars curr_succ_SUM
%%%%%%%%%%%%%%%%%%%%%%%%

% Extract the overlap between the image volume using sums by deleting the index less than the exponent
disp('Extract the overlap between the image volume using sums by deleting the index less than the exponent')
p_idx_mod(p_idx_mod(:,1)<exponent,:) = [];

% Create a variable curr_succ_variable where the columns have been annotated as
% (a) Column 1: Current Indices 
% (b) Column 2: Successive Indices 
% (c) Column 3: Pixel overlap between the current and successive indices 
% (d) Column 4: The pixel count of Current Indices 
% (e) Column 5: The pixel count of Successive Indices 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (a) Column 1: Current Indices
curr_succ_overlap(:,1)    =  floor(p_idx_mod(:,1)/exponent);

% (b) Column 2: Successive Indices
curr_succ_overlap(:,2)    =  mod(p_idx_mod(:,1),exponent);

% (c) Column 3: Pixel overlap between the current and successive indices 
curr_succ_overlap(:,3)    =  p_idx_mod(:,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('relieve memory  ===>> p_idx_mod')
clearvars p_idx_mod
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Delete the zero values indices in each row 
disp('Delete the zero values indices in each row')
tic
curr_succ_overlap(curr_succ_overlap(:,1)==0,:)=[];
curr_succ_overlap(curr_succ_overlap(:,2)==0,:)=[];
toc

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the pixel indices and pixel count of the current sub-volume 
disp('Determine the pixel indices and pixel count of the current sub-volume')
tic
[p_idx_mod_1(:,2),p_idx_mod_1(:,1)] = hist(nonzeros(curr_vol(:)),unique(nonzeros(curr_vol(:))));
toc

% Determine the missing indices in the pixel indices count and include them in the list and update the table accordingly 
disp('Determine the missing indices in the pixel indices count and include them in the list and update the table accordingly')
tic
missing          =    find(~ismember(1:max(p_idx_mod_1(:,1)),p_idx_mod_1(:,1)));
missing          =    [missing' zeros(numel(missing'),1)];
p_idx_mod_1      =    vertcat(p_idx_mod_1, missing);
[~,sort_order]   =    sort(p_idx_mod_1(:,1));
p_idx_mod_1      =    p_idx_mod_1(sort_order,:);
toc

% (d) Column 4: The pixel count of Current Indices 
disp('Column 4: The pixel count of Current Indices')
curr_succ_overlap(:,4)    =   p_idx_mod_1([curr_succ_overlap(:,1)]',2); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('relieve memory  ===>>> p_idx_mod_1  &  curr_vol')
clearvars p_idx_mod_1   curr_vol
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clearvars p_idx_mod_2       missing

% Determine the pixel indices and pixel count of the successive sub-volume 
disp('Determine the pixel indices and pixel count of the successive sub-volume')
tic
[p_idx_mod_2(:,2),p_idx_mod_2(:,1)]  =  hist(nonzeros(succ_vol(:)),unique(nonzeros(succ_vol(:))));
toc

% Determine the missing indices in the pixel indices count and include them in the list and update the table accordingly 
disp('Determine the missing indices in the pixel indices count and include them in the list and update the table accordingly')
tic
missing          =    find(~ismember(1:max(p_idx_mod_2(:,1)),p_idx_mod_2(:,1)));
missing          =    [missing' zeros(numel(missing'),1)];
p_idx_mod_2      =    vertcat(p_idx_mod_2, missing);
[~,sort_order]   =    sort(p_idx_mod_2(:,1));
p_idx_mod_2      =    p_idx_mod_2(sort_order,:);
toc

% (e) Column 5: The pixel count of Current Indices 
disp('Column 5: The pixel count of Current Indices')
curr_succ_overlap(:,5)    = p_idx_mod_2([curr_succ_overlap(:,2)]',2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('relieve memory  ====>>> p_idx_mod_2  & succ_vol')
clearvars p_idx_mod_2     succ_vol 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[min_val,min_idx] =  min((curr_succ_overlap(:,4:5)),[],2);
curr_succ_overlap = [curr_succ_overlap min_idx  min_val];

% final table 
disp('final table')
table_fin = curr_succ_overlap(:,1:2);
table_fin(:,3) = 100*curr_succ_overlap(:,3)./curr_succ_overlap(:,7);
table_fin(:,4) = curr_succ_overlap(:,6);
table_fin(:,5) = curr_succ_overlap(:,3);

% Rearrange the fiber columns in the appropriate terms
% for the curr_succ_route SM/SS
disp('Rearrange the fiber columns in the appropriate terms')
unique_fib_idx = nonzeros(unique(table_fin(:,1)));
tic
MAIN_COUNT_probe12_TMP_1 = cell(1,numel(unique_fib_idx));
for ii = 1:numel(unique_fib_idx)
    sprintf('%s%d','Probe Arrange ',ii)
   MAIN_COUNT_probe12_TMP_1{ii} = (table_fin(ismember(table_fin(:,1), unique_fib_idx(ii)),:))';
end
toc

% for the curr_succ_route MS
disp('for the curr_succ_route MS')
unique_fib_idx = nonzeros(unique(table_fin(:,2)));
MAIN_COUNT_query12_TMP_1 = cell(1,numel(unique_fib_idx));
tic
for ii = 1:numel(unique_fib_idx)
   sprintf('%s%d','Query Arrange ',ii)
   MAIN_COUNT_query12_TMP_1{ii} = (table_fin(ismember(table_fin(:,2), unique_fib_idx(ii)),:))';
end
toc
MAIN_COUNT_query12_TMP_1(cell2mat(cellfun(@(x)size(x,2)==1,MAIN_COUNT_query12_TMP_1,'uni',0))) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Remove the supposed single matched pairs in the indexes

tic
disp('Extract the supposed single matched pairs in the indexes')
TMP_CELL      =   MAIN_COUNT_probe12_TMP_1(cell2mat(cellfun(@(x)size(x,2)==1,MAIN_COUNT_probe12_TMP_1,'uni',0)));
toc

tic
disp('Extract second row of the fibers')
TMP_DEL       =   cell2mat(cellfun(@(x)x(2,:),TMP_CELL,'uni',0));
toc

tic
disp('delete the sngle fibers from the probe variable')
MAIN_COUNT_probe12_TMP_1(cell2mat(cellfun(@(x)size(x,2)==1,MAIN_COUNT_probe12_TMP_1,'uni',0))) = [];
toc

tic
disp('Check if second row of the fibers is false')
MULTI_SUCC    =   cell2mat(cellfun(@(x)unique(x(2,:)),MAIN_COUNT_query12_TMP_1,'uni',0));
toc

tic
disp('Remove false single fibers')
TMP_CELL(ismember(TMP_DEL,MULTI_SUCC)) = [];
toc
% save the probe cell
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

varargout{              1}   =   MAIN_COUNT_probe12_TMP_1;
varargout{          end+1}   =   MAIN_COUNT_query12_TMP_1;
varargout{          end+1}   =   TMP_CELL;
end








