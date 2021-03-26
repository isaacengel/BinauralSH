function Z = mult3(X,Y)
% Page-wise multiplication of 3D array X (size = M x N x L) and matrix Y 
% (size = N x P). If pagemtimes is defined (for versions > R2020b) the 
% output is the same as when using that function.
% This is used, for instance, to multiply an HRTF (nfreqs x ndirs x nears)
% and a matrix of SH coefficients (ndirs x nsh).
%
% INPUT:
%   X = 3D array (M x N x L)
%   Y = 2D matrix (N x P)
%
% OUTPUT: 
%   Z = 3D array (M x P x L)
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% March 2021

if exist('pagemtimes','builtin')
    Z = pagemtimes(X,Y);
else
    Z = cat(3,X(:,:,1)*Y,X(:,:,2)*Y);
end