function Hnm = toSH_SpSub(H,N,az,el,w,Nmax,reg_eps)
% Transform HRTF to SH domain at order N using spatial subsampling [1].
% This is equivalent to virtual loudspeaker decoding (first introduced by 
% [2]) with "mode-matching" (pseudoinverse-based).
%
% SIMPLE USAGE EXAMPLE:
%   Hnm = toSH_SpSub(H,15,az,el);
%
% INPUT:
%   H = HRTF up to Nyquist frequency (nfreqs x ndirs x 2 ears)
%   N = target SH order
%   az = HRIR azimuth (ndirs x 1) in rad
%   el = HRIR elevation (ndirs x 1) in rad (0=top, pi/2=front)
%   w = quadrature weights (ndirs x 1); if empty, use pseudoinverse
%   Nmax = highest available HRTF order (def=35)
%   
% OUTPUT:
%   Hnm = HRTF's SH coefficients (nfreqs x (N+1)^2 x 2 ears)
%
% REFERENCES:
%   [1] Bernschütz, Benjamin, et al. "Binaural reproduction of plane waves
%       with reduced modal order." Acta Acustica united with Acustica 100.5
%       (2014): 972-983.
%   [2] McKeag, Adam, and David S. McGrath. "Sound field format to binaural
%       decoder with head tracking." Audio Engineering Society Convention
%       6r. Audio Engineering Society, 1996.
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021

%% Prepare parameters
if ~exist('Nmax','var') || isempty(Nmax)
    Nmax = 35;
end

%% Process
% First, get the SH-HRTF with highest available order
Ymax_inv = getYinv(Nmax,az,el,w,reg_eps);
Hnm_max = mult3(H,Ymax_inv);

% Then, resample to an appropriately sized sparse grid (Gauss)
sgrid = sofia_gauss(2*(N+1),N+1,0);
Ymax_s = AKsh(Nmax, [], sgrid(:,1)*180/pi, sgrid(:,2)*180/pi, 'real').';
Hs = mult3(Hnm_max,Ymax_s); % = cat(3,Hnm_max(:,:,1)*Ymax_s,Hnm_max(:,:,2)*Ymax_s);

% Finally, convert to SH again at the target low order
Ys_inv = getYinv(N,sgrid(:,1),sgrid(:,2));
Hnm = mult3(Hs,Ys_inv);

 