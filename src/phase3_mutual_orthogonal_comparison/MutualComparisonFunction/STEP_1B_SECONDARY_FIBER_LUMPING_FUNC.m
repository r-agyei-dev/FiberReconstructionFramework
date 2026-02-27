function [varargout] = STEP_1B_SECONDARY_FIBER_LUMPING_FUNC(varargin)

MAIN_COUNT_probe   =    varargin{1};
MAIN_COUNT_query   =    varargin{2};
MAIN_COUNT_single  =    varargin{3};


%% Probe to Query path 
MAIN_COUNT_probeTMP = MAIN_COUNT_probe;
MAIN_COUNT_queryTMP = MAIN_COUNT_query;
MAIN_COUNT_singleTMP = MAIN_COUNT_single;

% Sharing the rogue regions 
multi_probe_values = cell2mat(cellfun(@(x)unique(x(1,:)),MAIN_COUNT_probeTMP,'uni',0));
multi_query_values = cell2mat(cellfun(@(x)unique(x(2,:)),MAIN_COUNT_queryTMP,'uni',0));

%%
% MAIN_COUNT_probeTMP is a matrix where the first row is n-repeating number and the second row are distinct numbers  
% This suggests that the probing fiber is a large bodied region mapping smaller bodied region in the query volume 
% It is there expected that the 4th row of the variable be the "number 2" suggesting that all the fibers in 2 are significantly smaller 
% Consequence of 4th row variation:

% (a) If all the elements are 1 then this suggests a smaller volume fiber in the probe region intersecting n-number of larger fibers in the query volume 
% (b) If at least one of the fibers is 1 suggests that the fiber in question is mapped unto a larger fiber 

%%
% Similarly arguments follow for the MAIN_COUNT_queryTMP but in the converse

% PROBE TO QUERY:
for tt = 1:size(MAIN_COUNT_probeTMP,2)
sprintf('%s%d','PROBE TO QUERY run ',tt)
            % If the above statement is true then the fiber at location is

                % This suggests that query fibers are larger fibers sharing a smaller rogue probe fiber fiber           
         %% SCENARIO A: 
                % in this if all are 1 remove completely since these are small fibers sandwiched in large fibers 
                probe_val = unique(MAIN_COUNT_probeTMP{tt}(1,:));
                
                  if all(MAIN_COUNT_probeTMP{tt}(4,:)==1)            

                     % Remove all the probe columns associated with query values that are large latching fibers  
                     if any(ismember(MAIN_COUNT_probeTMP{tt}(2,:),multi_query_values))
                         
                         % Extract the terms that are in multi_query 
                         multi_query_ext1 = MAIN_COUNT_probeTMP{tt}(2,(ismember(MAIN_COUNT_probeTMP{tt}(2,:),multi_query_values)));
                         
                         % Find which of those terms are large latching fibers 
                         
                         multi_query_ext2 = zeros(1,numel(multi_query_ext1));
                         for mm = 1:numel(multi_query_ext1)                             
                             if all(MAIN_COUNT_queryTMP{cell2mat(cellfun(@(x)any(ismember(x(2,:),multi_query_ext1(mm))),MAIN_COUNT_queryTMP,'uni',0))}(4,:) == 1)   
                                multi_query_ext2(mm) = multi_query_ext1(mm);   
                             end    
                         end
                         
                         multi_query_ext2(multi_query_ext2==0) = [];
                         
                         % Remove from the probe_fiber mat 
                         if~isempty(multi_query_ext2)
                         MAIN_COUNT_probeTMP{tt}(:,ismember(MAIN_COUNT_probeTMP{tt}(2,:),multi_query_ext2)) = [];
                         end
                     end
                      
%                % Delete all the fibers that have probe columns that have the query rows in multiquery and are large                
%                MAIN_COUNT_probeTMP{tt} = {};     

          %%   SCENARIO B:
                  % In this scenario one all but some query fibers are lesser than probe fiber so the decision is made to keep the 
                  % that which satisfies the overlapping crietria and remove that which does not satisfies
                  
                 elseif any(MAIN_COUNT_probeTMP{tt}(4,:) == 1)

                              % (a) Keep probe that satisfies the larger fiber overlapping smaller volumes criteria
                              % (b) Delete probe term that does not satify the criteria
                              % (c) Delete associated terms in the query as well
                              
                              % Currently we shall keep the term that does not satisfy has a term in the query then delete from the probe else keep 
                              
                              % if term in probe that does not satisfy is not in query then keep else delete
                              query_find = MAIN_COUNT_probeTMP{tt}(2,(MAIN_COUNT_probeTMP{tt}(4,:)==1));
                              
                              for kk = 1:numel(query_find)
                                 if ismember(query_find(kk),multi_query_values)
                                    MAIN_COUNT_probeTMP{tt}(:,ismember(MAIN_COUNT_probeTMP{tt}(2,:),query_find(kk))) = [];
                                 end
                              end                             
                             
                  end
end

%% PROBE TO PROBE:
MAIN_COUNT_probeTMP(cell2mat(cellfun(@(x)isempty(x),MAIN_COUNT_probeTMP,'uni',0)))   =   [];
probe4 = [MAIN_COUNT_probeTMP{:}];
% check for competing probe fibers 
[idx_mod(:,2),idx_mod(:,1)] = hist(probe4(2,:),unique(probe4(2,:)));
probe_compete = idx_mod(idx_mod(:,2)>1,1);

if ~isempty(probe_compete)
for tt = 1:numel(probe_compete)
    
    tmp_idx = find(cell2mat(cellfun(@(x)any(ismember(x(2,:),probe_compete(tt))),MAIN_COUNT_probeTMP,'uni',0)));
%     tmp_mat = [MAIN_COUNT_probeTMP{tmp_idx}];
    tmp_mat = cell2mat(cellfun(@(x)x(:,ismember(x(2,:),probe_compete(tt))),MAIN_COUNT_probeTMP(tmp_idx),'uni',0));

   max_loc = find(tmp_mat(5,:) == max(tmp_mat(5,:)));
   
   if numel(max_loc) > 1
       del_loc = max_loc(1);
   else 
       del_loc = max_loc;
   end    
    
    
    for jj = 1:numel(tmp_idx)
        if MAIN_COUNT_probeTMP{tmp_idx(jj)}(5,ismember(MAIN_COUNT_probeTMP{tmp_idx(jj)}(2,:),probe_compete(tt)))~= max(tmp_mat(5,:)) || ...
          (MAIN_COUNT_probeTMP{tmp_idx(jj)}(5,ismember(MAIN_COUNT_probeTMP{tmp_idx(jj)}(2,:),probe_compete(tt)))== max(tmp_mat(5,:)) && del_loc ~= jj)
      
          MAIN_COUNT_probeTMP{tmp_idx(jj)}(:,ismember(MAIN_COUNT_probeTMP{tmp_idx(jj)}(2,:),probe_compete(tt))) = [];
        end
    end
%     
end
end

%% QUERY TO PROBE:
MAIN_COUNT_probeTMP(cell2mat(cellfun(@(x)isempty(x),MAIN_COUNT_probeTMP,'uni',0)))   =   [];
MAIN_COUNT_queryTMP(cell2mat(cellfun(@(x)isempty(x),MAIN_COUNT_queryTMP,'uni',0)))   =   [];

multi_query_values = cell2mat(cellfun(@(x)unique(x(2,:)),MAIN_COUNT_queryTMP,'uni',0));
multi_probe_values = cell2mat(cellfun(@(x)unique(x(1,:)),MAIN_COUNT_probeTMP,'uni',0));

% QUERY TO PROBE:
for tt = 1:size(MAIN_COUNT_queryTMP,2)
sprintf('%s%d','QUERY TO PROBE run ',tt)
            % If the above statement is true then the fiber at location is

                % This suggests that query fibers are larger fibers sharing a smaller rogue probe fiber fiber so
                % (a) Delete the current probe cell and keep the combination in the better of the query cells sharing that rogue region
                % tt = 5;
           
         %% SCENARIO A: 
                % In this scenario smaller rogue fibers map to larger query fibers and a decision to allocate to which query faction is made 
                % Factions may be multi-linked or may be singly-linked
                % A1 yields: Probe_cells are deleted and Query cells are updated to maintain one that best matches probe fiber idx 
                % A2 yields: Probe_cells are maintained to keep the singly-linked and all query cells associated with the multi link are updated 
                query_val = unique(MAIN_COUNT_queryTMP{tt}(2,:));
                
                  if all(MAIN_COUNT_queryTMP{tt}(4,:)==2)
                      
                 % Remove all the query columns associated with probe values that are large latching fibers  
                     if any(ismember(MAIN_COUNT_queryTMP{tt}(1,:),multi_probe_values))
                         
                         % Extract the terms that are in multi_probe 
                         multi_probe_ext1 = MAIN_COUNT_queryTMP{tt}(1,(ismember(MAIN_COUNT_queryTMP{tt}(1,:),multi_probe_values)));
                         
                         % Find which of those terms are large latching fibers 
                         
                         multi_probe_ext2 = zeros(1,numel(multi_probe_ext1));
                         for mm = 1:numel(multi_probe_ext1)                             
                             if all(MAIN_COUNT_probeTMP{cell2mat(cellfun(@(x)any(ismember(x(1,:),multi_probe_ext1(mm))),MAIN_COUNT_probeTMP,'uni',0))}(4,:) == 2)   
                                multi_probe_ext2(mm) = multi_probe_ext1(mm);   
                             end    
                         end
                         
                         multi_probe_ext2(multi_probe_ext2==0) = [];
                         
                         % Remove from the query_fiber mat 
                         if~isempty(multi_probe_ext2)
                         MAIN_COUNT_queryTMP{tt}(:,ismember(MAIN_COUNT_queryTMP{tt}(1,:),multi_probe_ext2)) = [];
                         end
                     end    
%                  % Delete the query fiber cell 
%                  MAIN_COUNT_queryTMP{tt} = {};   
           
          %%   SCENARIO B:
                  % In this scenario one all but some query fibers are lesser than probe fiber so the decision is made to keep the 
                  % that which satisfies the overlapping crietria and remove that which does not satisfies
                  
                elseif any(MAIN_COUNT_queryTMP{tt}(4,:) == 2)

                              % (a) Keep probe that satisfies the larger fiber overlapping smaller volumes criteria
                              % (b) Delete probe term that does not satify the criteria
                              % (c) Delete associated terms in the query as well
                              
                              % if term in probe that does not satisfy is not in query then keep else delete
                              probe_find = MAIN_COUNT_queryTMP{tt}(1,(MAIN_COUNT_queryTMP{tt}(4,:) == 2));
                              
                              for kk = 1:numel(probe_find)
                                 if ismember(probe_find(kk),multi_probe_values)
                                    MAIN_COUNT_queryTMP{tt}(:,ismember(MAIN_COUNT_queryTMP{tt}(1,:),probe_find(kk))) = [];
                                 end
                              end

                  end
end

MAIN_COUNT_probeTMP(cell2mat(cellfun(@(x)isempty(x),MAIN_COUNT_probeTMP,'uni',0)))   =   [];
MAIN_COUNT_queryTMP(cell2mat(cellfun(@(x)isempty(x),MAIN_COUNT_queryTMP,'uni',0)))   =   [];

%% QUERY TO QUERY
clearvars idx_mod 
MAIN_COUNT_queryTMP(cell2mat(cellfun(@(x)isempty(x),MAIN_COUNT_queryTMP,'uni',0)))   =   [];
query4 = [MAIN_COUNT_queryTMP{:}];
% check for competing probe fibers 
[idx_mod(:,2),idx_mod(:,1)] = hist(query4(1,:),unique(query4(1,:)));
query_compete = idx_mod(idx_mod(:,2)>1,1);

if ~isempty(query_compete)
for tt = 1:numel(query_compete)
    
    tmp_idx = find(cell2mat(cellfun(@(x)any(ismember(x(1,:),query_compete(tt))),MAIN_COUNT_queryTMP,'uni',0)));
%     tmp_mat = [MAIN_COUNT_queryTMP{tmp_idx}];
    tmp_mat = cell2mat(cellfun(@(x)x(:,ismember(x(1,:),query_compete(tt))),MAIN_COUNT_queryTMP(tmp_idx),'uni',0));
    
   max_loc = find(tmp_mat(5,:) == max(tmp_mat(5,:)));
   
   if numel(max_loc) > 1
       del_loc = max_loc(1);
   else 
       del_loc = max_loc;
   end
    
    % here one of them is equal 
    for jj = 1:numel(tmp_idx)
        
        if MAIN_COUNT_queryTMP{tmp_idx(jj)}(5,ismember(MAIN_COUNT_queryTMP{tmp_idx(jj)}(1,:),query_compete(tt)))~=max(tmp_mat(5,:)) || ...
           (MAIN_COUNT_queryTMP{tmp_idx(jj)}(5,ismember(MAIN_COUNT_queryTMP{tmp_idx(jj)}(1,:),query_compete(tt)))== max(tmp_mat(5,:)) && del_loc ~= jj)
       
           MAIN_COUNT_queryTMP{tmp_idx(jj)}(:,ismember(MAIN_COUNT_queryTMP{tmp_idx(jj)}(1,:),query_compete(tt))) = [];
                      
        end  
        
    end
       
end
end

MAIN_COUNT_probeTMP(cell2mat(cellfun(@(x)isempty(x),MAIN_COUNT_probeTMP,'uni',0)))   =   [];
MAIN_COUNT_queryTMP(cell2mat(cellfun(@(x)isempty(x),MAIN_COUNT_queryTMP,'uni',0)))   =   [];


varargout{        1}     =   MAIN_COUNT_probeTMP;
varargout{    end+1}     =   MAIN_COUNT_queryTMP;

%%
% M_COUNT_p = MAIN_COUNT_probeTMP;
% M_COUNT_q = MAIN_COUNT_queryTMP;

% %% FIND THE SHARED FIBERS BETWEEN THE PROBING AND THE QUERY FIBERS
% probe_cat        =   [MAIN_COUNT_probeTMP{:}];
% probe_check2_p   =   unique(probe_cat(1,:));
% query_cat        =   [MAIN_COUNT_queryTMP{:}];
% query_check2_p   =   unique(query_cat(1,:));
% 
% common_probe = probe_check2_p(ismember(probe_check2_p,query_check2_p));
% 
% for ww = 1:numel(common_probe)
%     
%     probeTMP_cell = MAIN_COUNT_probeTMP{cell2mat(cellfun(@(x)any(ismember(x(1,:),common_probe(ww))),MAIN_COUNT_probeTMP,'uni',0))};
%     queryTMP_cell_find = find(cell2mat(cellfun(@(x)any(ismember(x(1,:),common_probe(ww))),MAIN_COUNT_queryTMP,'uni',0)));
%     
%      if ~isempty(queryTMP_cell_find) && ~isempty(probeTMP_cell)
% 
%             % if the probe fiber is a one to many then delete all the columns of the query with the probe values (which of the probe will you delete!!!!)
%             if size(probeTMP_cell,2)>1
% 
%              MAIN_COUNT_queryTMP{queryTMP_cell_find}(:,ismember(MAIN_COUNT_queryTMP{queryTMP_cell_find}(1,:),common_probe(ww))) = [];
% 
%             elseif size(probeTMP_cell,2)==1
% 
%                        % here the overlapping criteria not satisfied
%                        % for satisfaction queryTMP(4,:) = 1 and probeTMP(4,:) = 2;
%                         if probeTMP_cell(4,:) == 1 || MAIN_COUNT_queryTMP{queryTMP_cell_find}(4,ismember(MAIN_COUNT_queryTMP{queryTMP_cell_find}(1,:),common_probe(ww))) == 2
% 
%                             % find which of the configuration is larger 
%                             if probeTMP_cell(4,:) == 2
%                                 large = 1;
%                                 large_percent = probeTMP_cell(3,:);
%                             elseif MAIN_COUNT_queryTMP{queryTMP_cell_find}(4,ismember(MAIN_COUNT_queryTMP{queryTMP_cell_find}(1,:),common_probe(ww))) == 1
%                                 large = 2;
%                                 large_percent = MAIN_COUNT_queryTMP{queryTMP_cell_find}(3,ismember(MAIN_COUNT_queryTMP{queryTMP_cell_find}(1,:),common_probe(ww)));
% 
%                             else 
%                                 large_percent = [];
%                             end
% 
%                             if ~isempty(large_percent)
%                                         % if the larger has a confluence of 30% of more then let the larger keep the association 
%                                         if large_percent >= 30 
% 
%                                                if large == 1
%                                                MAIN_COUNT_queryTMP{queryTMP_cell_find}(:,ismember(MAIN_COUNT_queryTMP{queryTMP_cell_find}(1,:),common_probe(ww))) = [];
%                                                elseif large == 2
%                                                MAIN_COUNT_probeTMP{cell2mat(cellfun(@(x)any(ismember(x(1,:),common_probe(ww))),MAIN_COUNT_probeTMP,'uni',0))} = [];
%                                                end
% 
%                                         % if the larger has a confluence of less than 30% then letthe most confluenced keep the association    
%                                         elseif large_percent < 30 
% 
%                                              if probeTMP_cell(5,:) >= MAIN_COUNT_queryTMP{queryTMP_cell_find}(5,ismember(MAIN_COUNT_queryTMP{queryTMP_cell_find}(1,:),common_probe(ww)))
%                                                 MAIN_COUNT_queryTMP{queryTMP_cell_find}(:,ismember(MAIN_COUNT_queryTMP{queryTMP_cell_find}(1,:),common_probe(ww))) = [];
%                                              else 
%                                                  MAIN_COUNT_probeTMP{cell2mat(cellfun(@(x)any(ismember(x(1,:),common_probe(ww))),MAIN_COUNT_probeTMP,'uni',0))} = [];
%                                              end
% 
%                                         end
%                             else 
% 
%                                lump_mat = [probeTMP_cell    MAIN_COUNT_queryTMP{queryTMP_cell_find}(:,ismember(MAIN_COUNT_queryTMP{queryTMP_cell_find}(1,:),common_probe(ww)))];
% 
%                                if numel(nonzeros(ismember(lump_mat(5,:),max(lump_mat(5,:))))) == 1
%                                   MAIN_COUNT_queryTMP{queryTMP_cell_find}(:,~ismember(MAIN_COUNT_queryTMP{queryTMP_cell_find}(5,:),max(lump_mat(5,:)))) = [];
%                                   probeTMP_cell(~ismember(probeTMP_cell(5,:),max(lump_mat(5,:)))) = [];
%                                     if isempty(probeTMP_cell)  
%                                   MAIN_COUNT_probeTMP{cell2mat(cellfun(@(x)any(ismember(x(1,:),common_probe(ww))),MAIN_COUNT_probeTMP,'uni',0))} = [];
%                                     end
%                                else
%                                   MAIN_COUNT_queryTMP{queryTMP_cell_find} = [];
%                                end
% 
% 
%                             end
% 
%                         % here the overlapping criteria is satisfied 
%                         else 
% 
%                                   % if the larger has a confluence of less than 30% then letthe most confluenced keep the association
%                                  if probeTMP_cell(5,:) >= MAIN_COUNT_queryTMP{queryTMP_cell_find}(5,ismember(MAIN_COUNT_queryTMP{queryTMP_cell_find}(1,:),common_probe(ww)))
%                                     MAIN_COUNT_queryTMP{queryTMP_cell_find}(:,ismember(MAIN_COUNT_queryTMP{queryTMP_cell_find}(1,:),common_probe(ww))) = [];
%                                  else 
%                                      MAIN_COUNT_probeTMP{cell2mat(cellfun(@(x)any(ismember(x(1,:),common_probe(ww))),MAIN_COUNT_probeTMP,'uni',0))} = [];
%                                  end
% 
%                         end
% 
%             end
% 
%         MAIN_COUNT_probeTMP(cell2mat(cellfun(@(x)isempty(x),MAIN_COUNT_probeTMP,'uni',0)))   =   [];
%         MAIN_COUNT_queryTMP(cell2mat(cellfun(@(x)isempty(x),MAIN_COUNT_queryTMP,'uni',0)))   =   [];   
%            
%      end
% end
% 
% %% FIND THE SHARED FIBERS BETWEEN THE QUERY FIBERS AND THE PROBING FIBERS
% probe_cat        =   [MAIN_COUNT_probeTMP{:}];
% probe_check2_q   =   unique(probe_cat(2,:));
% query_cat        =   [MAIN_COUNT_queryTMP{:}];
% query_check2_q   =   unique(query_cat(2,:));
% 
% common_query = query_check2_q(ismember(query_check2_q,probe_check2_q));
% 
% for ww = 1:numel(common_query)
%     
%     queryTMP_cell = MAIN_COUNT_queryTMP{cell2mat(cellfun(@(x)any(ismember(x(2,:),common_query(ww))),MAIN_COUNT_queryTMP,'uni',0))};
%     probeTMP_cell_find = find(cell2mat(cellfun(@(x)any(ismember(x(2,:),common_query(ww))),MAIN_COUNT_probeTMP,'uni',0)));
% 
%  if ~isempty(queryTMP_cell) && ~isempty(probeTMP_cell_find)  
%      
%                 % if the query fiber is a one to many then delete all the columns of the probe with the query values (which of the query will you delete!!!!)
%                 if size(queryTMP_cell,2)>1
% 
%                  MAIN_COUNT_probeTMP{probeTMP_cell_find}(:,ismember(MAIN_COUNT_probeTMP{probeTMP_cell_find}(2,:),common_query(ww))) = [];
% 
%                 end
% 
%                 if size(queryTMP_cell,2) == 1
% 
%                            % here the overlapping criteria not satisfied
%                            % for satisfaction queryTMP(4,:) = 1 and probeTMP(4,:) = 2;
%                             if queryTMP_cell(4,:) == 2 || MAIN_COUNT_probeTMP{probeTMP_cell_find}(4,ismember(MAIN_COUNT_probeTMP{probeTMP_cell_find}(2,:),common_query(ww))) == 1
% 
%                                 % find which of the configuration is larger 
%                                 if queryTMP_cell(4,:) == 1
%                                     large = 2;
%                                     large_percent = queryTMP_cell(3,:);
%                                 elseif MAIN_COUNT_probeTMP{probeTMP_cell_find}(4,ismember(MAIN_COUNT_probeTMP{probeTMP_cell_find}(2,:),common_query(ww))) == 2
%                                     large = 1;
%                                     large_percent = MAIN_COUNT_probeTMP{probeTMP_cell_find}(3,ismember(MAIN_COUNT_probeTMP{probeTMP_cell_find}(2,:),common_query(ww)));
%                                 else 
%                                     large_percent = [];
%                                 end
% 
%                                 if ~isempty(large_percent)
%                                             % if the larger has a confluence of 30% of more then let the larger keep the association 
%                                             if large_percent >= 30 
% 
%                                                    if large == 2
%                                                    MAIN_COUNT_probeTMP{probeTMP_cell_find}(:,ismember(MAIN_COUNT_probeTMP{probeTMP_cell_find}(2,:),common_query(ww))) = [];
%                                                    elseif large == 1
%                                                    MAIN_COUNT_queryTMP{cell2mat(cellfun(@(x)any(ismember(x(2,:),common_query(ww))),MAIN_COUNT_queryTMP,'uni',0))} = [];
%                                                    end
% 
%                                             % if the larger has a confluence of less than 30% then letthe most confluenced keep the association    
%                                             elseif large_percent < 30 
% 
%                                                  if queryTMP_cell(5,:) >= MAIN_COUNT_probeTMP{probeTMP_cell_find}(5,ismember(MAIN_COUNT_probeTMP{probeTMP_cell_find}(2,:),common_query(ww)))
%                                                     MAIN_COUNT_probeTMP{probeTMP_cell_find}(:,ismember(MAIN_COUNT_probeTMP{probeTMP_cell_find}(2,:),common_query(ww))) = [];
%                                                  else 
%                                                      MAIN_COUNT_queryTMP{cell2mat(cellfun(@(x)any(ismember(x(2,:),common_query(ww))),MAIN_COUNT_queryTMP,'uni',0))} = [];
%                                                  end
% 
% 
%                                             end
% 
%                                 else 
% 
%                                    lump_mat = [queryTMP_cell   MAIN_COUNT_probeTMP{probeTMP_cell_find}(:,ismember(MAIN_COUNT_probeTMP{probeTMP_cell_find}(2,:),common_query(ww)))];
% 
%                                    if numel(nonzeros(ismember(lump_mat(5,:),max(lump_mat(5,:))))) == 1
%                                       MAIN_COUNT_probeTMP{probeTMP_cell_find}(:,~ismember(MAIN_COUNT_probeTMP{probeTMP_cell_find}(5,:),max(lump_mat(5,:)))) = [];
%                                       queryTMP_cell(~ismember(queryTMP_cell(5,:),max(lump_mat(5,:)))) = [];
%                                         if isempty(queryTMP_cell)  
%                                       MAIN_COUNT_queryTMP{cell2mat(cellfun(@(x)any(ismember(x(2,:),common_query(ww))),MAIN_COUNT_queryTMP,'uni',0))} = [];
%                                         end
%                                    else
%                                       MAIN_COUNT_probeTMP{probeTMP_cell_find} = [];
%                                    end
% 
%                                 end
% 
%                             % here the overlapping criteria is satisfied 
%                             else 
% 
%                                       % if the larger has a confluence of less than 30% then letthe most confluenced keep the association
%                                      if queryTMP_cell(5,:) >= MAIN_COUNT_probeTMP{probeTMP_cell_find}(5,ismember(MAIN_COUNT_probeTMP{probeTMP_cell_find}(2,:),common_query(ww)))
%                                         MAIN_COUNT_probeTMP{probeTMP_cell_find}(:,ismember(MAIN_COUNT_probeTMP{probeTMP_cell_find}(2,:),common_query(ww))) = [];
%                                      else 
%                                          MAIN_COUNT_queryTMP{cell2mat(cellfun(@(x)any(ismember(x(2,:),common_query(ww))),MAIN_COUNT_queryTMP,'uni',0))} = [];
%                                      end
% 
%                             end
% 
%                 end
% 
%             MAIN_COUNT_queryTMP(cell2mat(cellfun(@(x)isempty(x),MAIN_COUNT_queryTMP,'uni',0)))   =   [];
%             MAIN_COUNT_probeTMP(cell2mat(cellfun(@(x)isempty(x),MAIN_COUNT_probeTMP,'uni',0)))   =   [];   
%  end         
% end
% 
% MAIN_COUNT_probeTMP(cell2mat(cellfun(@(x)isempty(x),MAIN_COUNT_probeTMP,'uni',0)))   =   [];
% MAIN_COUNT_queryTMP(cell2mat(cellfun(@(x)isempty(x),MAIN_COUNT_queryTMP,'uni',0)))   =   [];



