function [Hnm,fc] = toSH_FDTA(H,N,az,el,fs,w,fc,r,earAz,earEl)
% Transform HRTF to SH domain at order N using frequency-dependent 
% time-alignment [1]. Alignment is performed with the method in [2].
%
% SIMPLE USAGE EXAMPLE:
%   Hnm = toSH_FDTA(H,15,az,el,48000);
%
% INPUT:
%   H = HRTF up to Nyquist frequency (nfreqs x ndirs x 2 ears)
%   N = target SH order
%   az = HRIR azimuth (ndirs x 1) in rad
%   el = HRIR elevation (ndirs x 1) in rad (0=top, pi/2=front)
%   fs = sampling frequency in Hz
%   w = quadrature weights (ndirs x 1); if empty, use pseudoinverse
%   fc = cutoff frequency above which the HRIRs are aligned (def=2500)
%   r = head radius in m (def=0.085)
%   earAz = left/right ear azimuth (1 x 2) in rad (def=[pi/2, 3*pi/2])
%   earEl = left/right ear elevation (1 x 2) in rad (def = [pi/2, pi/2])
%   
% OUTPUT:
%   Hnm = HRTF's SH coefficients (nfreqs x (N+1)^2 x 2 ears)
%   fc = see above
%
% EXTERNAL DEPENDENCIES:
%   AKtools (www.ak.tu-berlin.de/aktools)
%
% REFERENCES:
%   [1] Zaunschirm, M., Schörkhuber, C., Höldrich, R., 2018. Binaural
%       rendering of Ambisonic signals by head-related impulse response
%       time alignment and a diffuseness constraint. The Journal of the
%       Acoustical Society of America 143, 3616–3627.
%       https://doi.org/10.1121/1.5040489
%   [2] Ben-Hur, Zamir, et al. "Efficient Representation and Sparse
%       Sampling of Head-Related Transfer Functions Using Phase-Correction
%       Based on Ear Alignment." IEEE/ACM Transactions on Audio, Speech,
%       and Language Processing 27.12 (2019): 2249-2262.
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021

%% Some parameters
if ~exist('fc','var') || isempty(fc)
    fc = 2500; % based on MagLS paper (Schorkhuber et al, 2018)
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
nfreqs = size(H,1);
f = linspace(0,fs/2,nfreqs); % frequency vector
c = 343; % speed of sound (m/s)
kr = 2*pi*f*r/c;

%% First, apply phase correction above cutoff frequency
p = earAlign(kr,az,el,earAz,earEl); % phase correction
fc_ind = find(f>fc,1,'first');
ind = fc_ind:nfreqs;
H(ind,:,:) = H(ind,:,:).*exp(-1i*p(ind,:,:)); % apply correction

%% Then, get the order-truncated SH-HRTF
Y = AKsh(N,[],az*180/pi,el*180/pi,'real').';
if ~exist('w','var') || ~isempty(w)
    Y_inv = 4*pi*w.*Y'; % if integrations weights are provided, use them
else
    Y_inv = pinv(Y); % if not, the pseudoinverse will do just fine
end
Hnm = mult3(H,Y_inv);

