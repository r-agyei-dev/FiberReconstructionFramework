function [Xsub,idx]=LIROWS(X,tol)
% LIROWS  Extracts a linearly independent subset of rows from matrix X.
%
% Purpose:
%   Identifies a numerically independent set of rows using QR decomposition
%   with column pivoting. This is useful before SVD or regression steps
%   where duplicate or dependent samples can cause instability.
%
% Inputs:
%   X   - Input matrix (rows are candidate samples)
%   tol - Rank tolerance used for numerical independence (default = 1e-10)
%
% Outputs:
%   Xsub - Matrix containing only the independent rows (transposed at end
%          to preserve legacy behavior of the original pipeline)
%   idx  - Indices of the selected independent rows relative to X
%
% Notes:
%   • If X is all zeros, no independent rows exist and empty outputs return.
%   • QR with pivoting is used for robust rank detection.
%   • The final transpose is intentional to maintain compatibility with
%     existing downstream code that expects this orientation.

     % Handle degenerate case: matrix has no nonzero entries
     if ~nnz(X) % X contains only zeros → no independent rows
         Xsub=[]; idx=[];
         return
     end

     % Set default tolerance if not provided
     if nargin<2, tol=1e-10; end

     % QR decomposition with column pivoting for rank detection
     [Q, R, E] = qr(X,0); 

     % Extract magnitude of R diagonal (handles vector edge case)
     if ~isvector(R)
         diagr = abs(diag(R));
     else
         diagr = R(1);   
     end

     % Numerical rank estimation
     r = find(diagr >= tol*diagr(1), 1, 'last');

     % Indices of independent rows (via pivoting order)
     idx = sort(E(1:r));

     % Extract the independent subset
     Xsub = X(:,idx);

     % Transpose maintained for compatibility with existing pipeline
     Xsub = Xsub';