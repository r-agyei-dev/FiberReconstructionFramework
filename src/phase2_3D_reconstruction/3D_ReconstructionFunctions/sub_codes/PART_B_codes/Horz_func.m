function [keep_idx_f] = Horz_func(keep_idx_f)

for ii = 1:size(keep_idx_f,2)   
    if ~isempty(keep_idx_f{ii}) &&  size(keep_idx_f{ii},1)>1       
        keep_idx_f{ii} = (keep_idx_f{ii})';        
    end
    
end 