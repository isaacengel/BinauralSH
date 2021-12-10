function Y_inv = getYinv(N,az,el,w,reg_eps)
% Get transformation matrix to obtain the SH coefficients of an HRTF.
% Available options:
%   1. Use provided integration weights (see w)
%   2. Use pseudoinverse method (least squares solution)
%   3. Use Tikhonov regularisation which is useful if there are missing
%       directions, e.g. at the bottom of the sphere [1]
%
% SIMPLE USAGE EXAMPLE:
%   Y_inv = getYinv(15,az,el);
%
% INPUT:
%   N = target SH order
%   az = HRIR azimuth (ndirs x 1) in rad
%   el = HRIR elevation (ndirs x 1) in rad (0=top, pi/2=front)
%   w = quadrature weights (ndirs x 1); if empty, use pseudoinverse
%   reg_eps = epsilon used for Tikhonov regularisation according to [1]; if
%       empty, set to 0 (don't regularise)
%
% OUTPUT:
%   Y_inv = transformation matrix (ndirs x (N+1)^2)
%
% REFERENCES:
%   [1] Duraiswaini, R., Dmitry N. Zotkin, and Nail A. Gumerov.
%       "Interpolation and range extrapolation of HRTFs." 2004 IEEE 
%       International Conference on Acoustics, Speech, and Signal
%       Processing. Vol. 4. IEEE, 2004.
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% August 2021

if ~exist('w','var')
    w = [];
end

if ~exist('reg_eps','var')
    reg_eps = 0; % if epsilon not provided, don't regularise
end

Y = getRealSHmatrix(N,az,el);

if ~isempty(w) && reg_eps == 0
    
    % If integration weights are provided, use them
%     Y_inv = 4*pi*w.*Y'; % not compatible with old Matlab versions
    Y_inv = mult2(4*pi*w,Y');
    
else
    
    if reg_eps == 0
        
        % If no regularisation requested, just apply the pseudoinverse
        Y_inv = pinv(Y); 
        
    else
        
        % If regualarisation requested, apply it according to [12]
        D = (1 + N*(N+1)) * eye((N+1)^2); % reg. matrix
        Y_inv = Y'*(Y*Y'+reg_eps*D)^(-1); % regularised weights
        
    end
end

