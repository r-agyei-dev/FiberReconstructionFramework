%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PRIMARY_WATERSHED_FUNCTION
% -------------------------------------------------------------------------
% PURPOSE:
% This function performs watershed segmentation on a 2D image and computes 
% the geometrical metadata of the segmented regions (ellipses).
%
% INPUTS:
%   Image_crop - Grayscale 2D image slice to be segmented
%   level      - Threshold for binarization (Otsu or manual, default ~0.5)
%   mask_value - Parameter controlling minima suppression for watershed
%
% OUTPUTS:
%   Stats - Structure containing segmented region metadata (PixelIdxList)
%
% DESCRIPTION:
%   1. Converts the image to binary using the specified threshold.
%   2. Computes the distance transform of the binary image.
%   3. Uses imextendedmin and imimposemin to control watershed minima.
%   4. Applies watershed segmentation to separate touching regions.
%   5. Removes small objects and extracts region properties.
% -------------------------------------------------------------------------

function [Stats] = PRIMARY_WATERSHED_FUNCTION(Image_crop,level,mask_value)

%% ------------------------------------------------------------------------
% Step 1: Binarize the image using the specified threshold
% -------------------------------------------------------------------------
Image_crop = imbinarize(Image_crop, level);  % Convert grayscale to binary

%% ------------------------------------------------------------------------
% Step 2: Compute distance transform for watershed
% -------------------------------------------------------------------------
D = -bwdist(~Image_crop, 'euclidean');      % Distance transform (negative for watershed minima)

% figure, imshow(D, []);                    % Uncomment to visualize distance transform

%% ------------------------------------------------------------------------
% Step 3: Impose minima and apply watershed
% -------------------------------------------------------------------------
Ld2 = watershed(imimposemin(D, imextendedmin(D, mask_value)));  % Watershed segmentation
% imextendedmin identifies regional minima below mask_value threshold
% imimposemin forces these minima onto the distance transform

Image_crop(Ld2 == 0) = 0;                    % Set watershed lines to background

%% ------------------------------------------------------------------------
% Step 4: Extract properties of segmented regions
% -------------------------------------------------------------------------
Stats = regionprops(bwareaopen(Image_crop, 5), 'PixelIdxList');  
% bwareaopen removes small objects (<5 pixels)
% PixelIdxList contains linear indices of pixels in each region

end