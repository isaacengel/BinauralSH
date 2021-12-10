function [Hnm,fc] = toSH_SpSubMod(H,N,az,el,fs,w,fc,r,earAz,earEl,Nmax,reg_eps)
% Combination of spatial subsampling [1] and frequency-dependent time 
% alignment [2], as suggested in [3].
%
% SIMPLE USAGE EXAMPLE:
%   Hnm = toSH_SpSubMod(H,15,az,el,48000);
%
% INPUT:
%   H = HRTF up to Nyquist frequency (nfreqs x ndirs x 2 ears)
%   N = target SH order
%   az = HRIR azimuth (ndirs x 1) in rad
%   el = HRIR elevation (ndirs x 1) in rad (0=top, pi/2=front)
%   fs = sampling frequency in Hz
%   w = quadrature weights (ndirs x 1); if empty, use pseudoinverse
%   fc = cutoff frequency above which the HRIRs are aligned (def=1500)
%   r = head radius in m (def=0.085)
%   earAz = left/right ear azimuth (1 x 2) in rad (def=[pi/2, 3*pi/2])
%   earEl = left/right ear elevation (1 x 2) in rad (def = [pi/2, pi/2])
%   Nmax = highest available HRTF order (def=35)
%   
% OUTPUT:
%   Hnm = HRTF's SH coefficients (nfreqs x (N+1)^2 x 2 ears)
%   fc = see above
%
% REFERENCES:
%   [1] Zaunschirm, M., Schörkhuber, C., Höldrich, R., 2018. Binaural
%       rendering of Ambisonic signals by head-related impulse response
%       time alignment and a diffuseness constraint. The Journal of the
%       Acoustical Society of America 143, 3616–3627.
%       https://doi.org/10.1121/1.5040489
%   [2] Bernschütz, Benjamin, et al. "Binaural reproduction of plane waves
%       with reduced modal order." Acta Acustica united with Acustica 100.5
%       (2014): 972-983.
%   [3] McKenzie, T., Murphy, D., Kearney, G., 2019. An evaluation of
%       pre-processing techniques for virtual loudspeaker binaural
%       ambisonic rendering, in: EAA Spatial Audio Signal Processing
%       Symposium. Paris, France, pp. 149–154.
%       https://doi.org/10.25836/sasp.2019.09
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021

%% Some parameters
if ~exist('w','var')
    w = [];
end
if ~exist('fc','var') || isempty(fc)
    fc = 2500; % same value used in [1]
end
if ~exist('r','var') || isempty(r)
    r = 0.085;
end
if ~exist('earAz','var') || isempty(earAz)
    earAz = [pi/2, 3*pi/2];
end
if ~exist('earEl','var') || isempty(earEl)
    earEl = [pi/2, pi/2];
end
if ~exist('Nmax','var') || isempty(Nmax)
    Nmax = 35;
end
nfreqs = size(H,1);
f = linspace(0,fs/2,nfreqs); % frequency vector
c = 343; % speed of sound (m/s)
kr = 2*pi*f*r/c;

%% First, apply phase correction above cutoff frequency
p = earAlign(kr,az,el,earAz,earEl); % phase correction
fc_ind = find(f>fc,1,'first');
ind = fc_ind:nfreqs;
H(ind,:,:) = H(ind,:,:).*exp(-1i*p(ind,:,:)); % apply correction

%% Then, apply spatial subsampling
Hnm = toSH_SpSub(H,N,az,el,w,Nmax,reg_eps);
