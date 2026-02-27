%                                                                                         Function was created by Ronald F Agyei.
%                                                                                         Date: 10/04/2018

% This function was created to efficiently manages the .mat files from loadeded images 
% Note that Large .mat files are memory intensive

function [varargout] = Tiff_Fast_Loader(varargin)

fileFolder = varargin{1};
type_ext   = varargin{2};

cd(fileFolder);
filePattern = fullfile(fileFolder,type_ext);
dirOutput = dir(filePattern);
fileNames = {dirOutput.name}';
curr_cell = cell(1,numel(fileNames));

tic 
for jj = 1:numel(fileNames)
curr_cell{jj} = imread(fileNames{(jj)});        
end
toc

tic 
curr_cell = vertcat(curr_cell{:});
[r,c] = size(curr_cell);
nlay  = numel(fileNames);
curr_cell   = permute(reshape(curr_cell',[c,r/nlay,nlay]),[2,1,3]);
toc 

varargout{       1}  =  curr_cell;
end