function [Hnm,fc] = toSH_MagSun(H,N,az,el,fs,w,fc,k,r,reg_eps)
% Transform HRTF to SH domain at order N with the Magnitude-optimised
% reconstruction method from Sun's thesis [1].
%
% SIMPLE USAGE EXAMPLE:
%   Hnm = toSH_MagSun(H,15,az,el,48000);
%
% INPUT:
%   H = HRTF up to Nyquist frequency (nfreqs x ndirs x 2 ears)
%   N = target SH order
%   az = HRIR azimuth (ndirs x 1) in rad
%   el = HRIR elevation (ndirs x 1) in rad (0=top, pi/2=front)
%   fs = sampling frequency in Hz
%   w = quadrature weights (ndirs x 1); if empty, use pseudoinverse
%   fc = cutoff frequency in Hz above which phase is "disregarded" in
%       favour of magnitude; if empty (default), use aliasing frequency
%   k = length of the transition band in octaves (def=1). E.g.
%       if fc=623.89 Hz, smooth between 349.65 and 882.31 Hz. If k==0, 
%       don't smooth.
%   r = head radius in m (def=0.0875)
%   
% OUTPUT:
%   Hnm = HRTF's SH coefficients (nfreqs x (N+1)^2 x 2 ears)
%   fc = see above
%
% REFERENCES:
%   [1] Sun, David. "Generation and perception of three-dimensional sound
%       fields using Higher Order Ambisonics." (2013).
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% September 2021

%% Some parameters
if ~exist('r','var') || isempty(r)
    r = 0.0875; % default head radius
end
if ~exist('fc','var') || isempty(fc)
    c = 343; % speed of sound (m/s)
    fc = N*c/(2*pi*r); % if fc not provided, use aliasing frequency
end
if ~exist('k','var') || isempty(k)
    k = 1; % default smoothing = one octave (half to each side)
end
nfreqs = size(H,1);
f = linspace(0,fs/2,nfreqs).'; % frequency vector

%% Get the SH matrices

Hmag = abs(H);
Hphase = unwrap(angle(H));
HphaseAvg = mean(Hphase,2); %HphaseAvg = mult3(Hphase,ones(ndirs))/ndirs;
alpha = 0.5 + k/log(2) * log(f/fc); % 0 for fc1, 1 for fc2
alpha = min(max(alpha,0),1); % keep between 0 and 1
HphaseMod = alpha.*HphaseAvg + (1-alpha).*Hphase;
Htarg = Hmag.*exp(1i*HphaseMod);
Y_inv = getYinv(N,az,el,w,reg_eps);
Hnm = mult3(Htarg,Y_inv); 
