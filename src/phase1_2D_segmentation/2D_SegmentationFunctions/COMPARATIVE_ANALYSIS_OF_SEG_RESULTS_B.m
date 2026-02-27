%--------------------------------------------------------------------------
% COMPARATIVE_ANALYSIS_OF_SEG_RESULTS_B
%
% PURPOSE:
%   Compares segmentation results obtained from an iteratively sharpened
%   image with those from the original (OTSU-thresholded) image. The
%   function restores missing regions and replaces certain one-to-one
%   matches with the original configurations when appropriate.
%
% INPUTS:
%   Stats_Optimized -> regionprops metadata from sharpened segmentation
%   Stats_Crude     -> regionprops metadata from OTSU segmentation
%   size_length     -> size vector of the grayscale image
%
% OUTPUT:
%   varargout{1}    -> updated Stats structure after reconciliation
%--------------------------------------------------------------------------

function [varargout] = COMPARATIVE_ANALYSIS_OF_SEG_RESULTS_B(varargin)

%--------------------------------------------------------------------------
% Unpack inputs
%--------------------------------------------------------------------------
Stats_Optimized       = varargin{1};
Stats_Crude           = varargin{2};
size_length           = varargin{3};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 1: Build labeled maps for sharpened and crude segmentations
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% (a) Extract pixel index lists from regionprops outputs
Pix_idx_curr   = {Stats_Optimized.Pixel_IDX_List};
Pix_idx_consec = {Stats_Crude.Pixel_IDX_List};

% (b) Initialize label images and duplicate trackers
Temp_sharp_curr   = zeros(size_length);   % label map (sharpened)
Temp_sharp_consec = zeros(size_length);   % label map (crude)
double_ganger_curr = zeros(1,length(Pix_idx_curr));   % duplicate flags (curr)
double_ganger_succ = zeros(1,length(Pix_idx_consec)); % duplicate flags (consec)

% Populate sharpened label image
for tt = 1:length(Pix_idx_curr)
    if Temp_sharp_curr(Pix_idx_curr{tt}) == 0
        Temp_sharp_curr(Pix_idx_curr{tt}) = tt;
    else
        % mark regions that collide (duplicate ownership)
        double_ganger_curr(tt) = tt;
    end
end

% Populate crude label image
for tt = 1:length(Pix_idx_consec)
    if Temp_sharp_consec(Pix_idx_consec{tt}) == 0
        Temp_sharp_consec(Pix_idx_consec{tt}) = tt;
    else
        % mark regions that collide (duplicate ownership)
        double_ganger_succ(tt) = tt;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% STEP 2: Determine region correspondences and update structures
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Run matching analysis between sharpened and crude label maps
[Curr_Cell,~,Consec_uncaptured_cell,~] = ...
    EPAD_FUNCTION_CODE_MASTER_VERSION_2(Temp_sharp_curr,Temp_sharp_consec);

% Extract one-to-one matches (single correspondences)
Curr_Single_Cell = Curr_Cell( ...
    cell2mat(cellfun(@(x) size(x,2)==1, Curr_Cell,'uni',0)));

% Separate indices for deletion (from optimized) and addition (from crude)
Curr_single_delete = cell2mat(cellfun(@(x)x(1),Curr_Single_Cell,'uni',0));
Succ_single_add    = cell2mat(cellfun(@(x)x(2),Curr_Single_Cell,'uni',0));

% Remove selected optimized regions and known duplicates
Stats_Optimized([Curr_single_delete ...
    double_ganger_curr(double_ganger_curr~=0)]) = [];

% Remove crude regions that were duplicates
Consec_uncaptured_cell( ...
    ismember(Consec_uncaptured_cell, ...
    double_ganger_succ(double_ganger_succ~=0))) = [];

% Append recovered crude regions to the optimized structure
Stats_Optimized = [Stats_Optimized; ...
                   Stats_Crude(Succ_single_add); ...
                   Stats_Crude(Consec_uncaptured_cell)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Output
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
varargout{1} = Stats_Optimized;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%