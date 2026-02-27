
function [varargout] = MUTUAL_COMPARISON_CODE_ANALYZERS(varargin)

ZXY_YZX_SM        =   varargin{1};
YZX_ZXY_SM        =   varargin{2};
ZXY_YZX_SS        =   varargin{3};


% Extraction of the fibers from the comparison 
Curr_Idx     =   cell2mat(cellfun(@(x)unique(x(1,:)),ZXY_YZX_SM,'uni',0));
Succ_Idx     =   cell2mat(cellfun(@(x)unique(x(2,:)),YZX_ZXY_SM,'uni',0));
single_Idx   =   cell2mat(cellfun(@(x)[x(x(4)) ; x(4)],ZXY_YZX_SS,'uni',0));

Curr_Idx = unique([Curr_Idx  single_Idx(1,single_Idx(2,:) == 1)]);
Succ_Idx = unique([Succ_Idx  single_Idx(1,single_Idx(2,:) == 2)]);

varargout{         1}  =   Curr_Idx;
varargout{     end+1}  =   Succ_Idx;
end





