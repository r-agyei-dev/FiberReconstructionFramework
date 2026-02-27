

function [cell_var] = Extract_LIDX_fun(VOL)


divisor = 100;
slice_depth = numel(VOL);

clearvars VOL

slice = single(1:round((slice_depth/divisor)):slice_depth);

cell_var = cell(1,numel(slice));
 
 for uu = 1:numel(slice)
     if uu == 1
         slice_tmp = 1:slice(uu+1)-1;
     elseif uu > 1 && uu < numel(slice)
         slice_tmp = slice(uu):slice(uu+1)-1;
     elseif uu == numel(slice)
         slice_tmp = slice(uu):slice_depth;
     end  
     cell_var{uu} = slice_tmp;
 end

end