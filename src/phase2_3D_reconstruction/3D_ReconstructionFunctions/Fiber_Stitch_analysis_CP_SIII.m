% Instead of saving and reloading intermediate files, this function
% performs the EPAD-based stitching analysis directly in memory and
% returns the necessary outputs.

function [varargout] = Fiber_Stitch_analysis_CP_SIII(varargin)

% =========================================================================
% INPUT PARSING
% =========================================================================
Linear_Index_Final           =     varargin{1}; % Final linear indices per fiber
height_variable_MAIN         =     varargin{2}; % Slice index per pixel for each fiber
Pixels_NON_BIN               =     varargin{3}; % Mask of non-connected pixels
size_length                  =     varargin{4}; % Size of the 3D volume

%%
tic

        disp('Begin Image instantiations')

        % -----------------------------------------------------------------
        % Allocate working volumes:
        % IMAGE_VOL_1 → fiber ID labels
        % IMAGE_VOL_2 → slice index labels
        % -----------------------------------------------------------------

        IMAGE_VOL_1 = zeros(size_length);
        IMAGE_VOL_2 = zeros(size_length);        

        % Preallocate containers for EPAD processing
        curr_values = cell(1,length(Linear_Index_Final));
        succ_values = cell(1,length(Linear_Index_Final));
        
        % -----------------------------------------------------------------
        % Populate volumes with fiber labels and slice information
        % -----------------------------------------------------------------
        for ii = 1:length(Linear_Index_Final)

%             sprintf('%s%d','image instiations ',ii)

            % -------------------------------------------------------------
            % Label current fiber ID into IMAGE_VOL_1
            % -------------------------------------------------------------
            IMAGE_VOL_1(Linear_Index_Final{ii}) = ii;                    
            curr_values{ii} = IMAGE_VOL_1(Linear_Index_Final{ii});
            
            % -------------------------------------------------------------
            % Label slice index into IMAGE_VOL_2
            % -------------------------------------------------------------
            IMAGE_VOL_2(Linear_Index_Final{ii}) = height_variable_MAIN{ii};
            
            % -------------------------------------------------------------
            % Remove non-connected regions by zeroing them out
            % -------------------------------------------------------------
            % disp('Set the regions with the non_connected regions to zero')
            IMAGE_VOL_2(Linear_Index_Final{ii}(Pixels_NON_BIN(Linear_Index_Final{ii}))) = 0;
            succ_values{ii} = IMAGE_VOL_2(Linear_Index_Final{ii});

        end

        disp('End Image instantiations')
        
        % -----------------------------------------------------------------
        % Memory cleanup (intentionally preserved)
        % -----------------------------------------------------------------
        clearvars height_variable_MAIN Linear_Index_Final Pixels_NON_BIN
        
        % Additional cleanup (kept for compatibility)
        clearvars Pixels_NON_BIN
        
        % -----------------------------------------------------------------
        % Execute EPAD algorithm to resolve overlapping slices
        % -----------------------------------------------------------------
        disp('Use the EPAD Algorithm and solve')
        
        % Construct exponent used internally by EPAD encoding
        exponent = 10^(numel(num2str(max(IMAGE_VOL_2(:)))) + 2);

        % Run EPAD fiber stitching analysis
        [curr_succ_multiple,split_array_consec_int,split_array_slice] = ...
            EPAD_FUNCTION_CODE_Fiber_Stitch_ver(curr_values,succ_values,exponent);

        % Free large temporary volumes
        clearvars IMAGE_VOL_1 IMAGE_VOL_2     
        
toc

% =========================================================================
% OUTPUTS
% =========================================================================
varargout{           1}   =   split_array_slice;        % Slice-based split info
varargout{       end+1}   =   split_array_consec_int;   % Consecutive intersection info
end