% =========================================================================
% SHARP_SEGMENTOR_ADJUSTOR
% =========================================================================
% This function performs image sharpening, intensity adjustment, and
% iterative watershed segmentation on a raw grayscale image slice.
%
% INPUTS:
%   Image_sharpen      - Grayscale image slice to process
%   thresh_val_curr    - Threshold value for intensity adjustment
%
% OUTPUTS:
%   Stats_UNSHARP      - Structure array containing ellipse metadata for
%                        the segmented regions
%   Image_sharpen      - Adjusted and sharpened image
% =========================================================================

function [Stats_UNSHARP, Image_sharpen] = SHARP_SEGMENTOR_ADJUSTOR(Image_sharpen, thresh_val_curr)    

%% ------------------------------------------------------------------------
% Step 1: Image normalization and intensity adjustment
% -------------------------------------------------------------------------
% Normalize image to [0,1] and apply contrast adjustment based on mean and
% standard deviation using the user-specified threshold
Image_sharpen = Image_sharpen / 65535;
avg = mean2(Image_sharpen); 
sigma = std2(Image_sharpen);
Image_sharpen = uint16(65535 * imadjust(Image_sharpen, ...
    [avg - thresh_val_curr * sigma, avg + thresh_val_curr * sigma], []));
figure, imshow(Image_sharpen);

%% ------------------------------------------------------------------------
% Step 2: Image cleaning (optional edge removal or cropping)
% -------------------------------------------------------------------------
% Currently, no cleaning is applied. Placeholder for future preprocessing
Clean_Image_UpdateS = Image_sharpen;

%% ------------------------------------------------------------------------
% Step 3: Initial watershed segmentation
% -------------------------------------------------------------------------
Image_crop = Clean_Image_UpdateS;
fib_diam_thresh = 10.38;        % Threshold for minor axis length
level_orig = 0.5;               % Threshold for binarization
mask_value_orig = 1.1;          % Initial mask value for watershed

[~, Stats] = WATERSHED_FUNC_UPDATE_SMART(Image_crop, level_orig, mask_value_orig);

%% ------------------------------------------------------------------------
% Step 4: Iterative watershed segmentation for oversized regions
% -------------------------------------------------------------------------
if ~isempty(Stats)
    
    % Define descending mask values for iteration
    mask_value_iterate = sort(0.2:0.05:1.0, 'descend');
    
    % Preallocate Stats storage
    Stats_cell{1, numel(mask_value_iterate)} = [];
    Stats_cell{1} = Stats;
    Stats_to_keep{1, numel(mask_value_iterate)} = [];
    
    % Iterate over mask values
    for jj = 1:numel(mask_value_iterate)
        Minor_Vec = [Stats_cell{jj}.MinorAxisLength];
        if ~isempty(Minor_Vec > fib_diam_thresh)   % Check for regions above threshold
            [Stats_to_keep{jj}, Stats_cell{jj+1}, ~] = ITERATE_WATERSHED_UPDATE_SMART( ...
                Image_crop, level_orig, mask_value_iterate(jj), Stats_cell{jj});
        else
            break
        end
    end
    
    % Lump all intermediate results into one array
    Stats_to_keep = STATS_LUMP_UPDATE(Stats_to_keep);
    Stats_cell_lumped = STATS_LUMP_UPDATE(Stats_cell);
    
    % Merge kept stats with last iteration
    if isempty(Stats_to_keep)
        Stats_Update_OUT = Stats_cell_lumped;
    else
        Stats_Update_OUT = [Stats_to_keep; Stats_cell{end}];
    end
    
    Stats = Stats_Update_OUT;
else
    Stats = [];
end

%% ------------------------------------------------------------------------
% Step 5: Recompute ellipse statistics for the segmented regions
% -------------------------------------------------------------------------
if ~isempty(Stats)
    Stats_FINAL_UN_Sharp = FINAL_STATS_FORMULATOR_SMART(Stats, Image_crop);
    Stats_UNSHARP = Stats_FINAL_UN_Sharp;
else
    Stats_UNSHARP = [];
end

% Optionally, save results for later use
% save('Stats_UNSHARP.mat', 'Stats_FINAL_UN_Sharp', '-v7.3')
