% ITERATIVE_IM_SHARP_WATERSHED_SEG_FUNCTION
% -------------------------------------------------------------------------
% PURPOSE:
% This function performs iterative image preprocessing (adjustment and sharpening)
% followed by watershed segmentation on 2D image slices. It outputs both
% the sharpened image and the segmentation metadata.
%
% INPUTS:
%   Image_sharpen  - Grayscale 2D image slice
%   d_amt_rad_adj  - Array of parameters for imadjust, imsharpen, and watershed
%                    [Amount, Radius, Threshold, Imadjust_max, Imadjust_min, Watershed_level]
%   Adjust         - Flag indicating whether image adjustment/sharpening is needed (1 = yes)
%   size_length_2D - Size of the 2D image slice ([rows, columns])
%   columnsInImage - Column meshgrid (for pixel indexing)
%   rowsInImage    - Row meshgrid (for pixel indexing)
%
% OUTPUTS:
%   Stats_WSEG_Reg - Structure containing statistics of watershed-segmented regions
%   Image_sharpen  - Preprocessed (adjusted/sharpened) image
% -------------------------------------------------------------------------

function [Stats_WSEG_Reg,Image_sharpen] = ITERATIVE_IM_SHARP_WATERSHED_SEG_FUNCTION(varargin)

%% ------------------------------------------------------------------------
% Extract inputs
% -------------------------------------------------------------------------
Image_sharpen    = varargin{1};
d_amt_rad_adj    = varargin{2};
Adjust           = varargin{3};
size_length_2D   = varargin{4};
columnsInImage   = varargin{5};
rowsInImage      = varargin{6};

%% ------------------------------------------------------------------------
% Image adjustment and sharpening (optional)
% -------------------------------------------------------------------------
if Adjust == 1
    % Adjust image intensity range
    Image_sharpen = imadjust(Image_sharpen, [d_amt_rad_adj(5) d_amt_rad_adj(4)], []);
    % Sharpen image using specified radius, amount, and threshold
    Image_sharpen = imsharpen(Image_sharpen, 'Radius', d_amt_rad_adj(2), ...
                                           'Amount', d_amt_rad_adj(1), ...
                                           'Threshold', d_amt_rad_adj(3));
end

% figure, imshow(Image_sharpen) % Uncomment to visualize preprocessed image

%% ------------------------------------------------------------------------
% Perform watershed segmentation
% -------------------------------------------------------------------------
mask_value_orig = 1.1;

% PRIMARY_WATERSHED_FUNCTION outputs ellipse statistics of segmented regions
[Stats] = PRIMARY_WATERSHED_FUNCTION(uint16(Image_sharpen), d_amt_rad_adj(6), mask_value_orig);

%% ------------------------------------------------------------------------
% Recompute ellipse-based statistics for segmented regions
% -------------------------------------------------------------------------
if ~isempty(Stats)
    % Convert pixel-based regions to ellipse-based metrics (area, orientation, centroid, etc.)
    [Stats_WSEG_Reg] = FINAL_STATS_FORMULATOR(Stats, size_length_2D, columnsInImage, rowsInImage);
else
    Stats_WSEG_Reg = [];
end

end