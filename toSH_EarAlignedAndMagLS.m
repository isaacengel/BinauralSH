function [Hnm,fc] = toSH_EarAlignedAndMagLS(H,N,az,el,fs,w,fc,frac,r,earAz,earEl)
% Transform HRTF to SH domain at order N by using ear aligning first [1]
% and frequency-dependent optimisation above a certain cutoff frequency
% (MagLS [2,3]). The idea is to align the HRIRs first to reduce the
% effective SH order, so we don't have to worry about linear phase and
% ITDs, and then transform to SH domain, but optimising for magnitude above
% a certain frequency (e.g. duplex theory).
%
% Important: the method works better if the HRIRs are aligned at t=0. This
% is taken care of in the function 'toSH'.
%
% To undo the alignment after SH interpolation, do:
%
%   H_aligned = Hnm * Y; % Y = SH coeffs for desired direction
%   H = H_aligned.*exp(1i*p); % p = second output of this function
%   h = iffth(H);
%
% SIMPLE USAGE EXAMPLE:
%   Hnm = toSH_EarAlignedAndMagLS(H,15,az,el,48000);
%
% INPUT:
%   H = HRTF up to Nyquist frequency (nfreqs x ndirs x 2 ears)
%   N = target SH order
%   az = HRIR azimuth (ndirs x 1) in rad
%   el = HRIR elevation (ndirs x 1) in rad (0=top, pi/2=front)
%   fs = sampling frequency in Hz
%   w = quadrature weights (ndirs x 1); if empty, use pseudoinverse
%   fc = cutoff frequency in Hz above which phase is "disregarded" in
%   favour of magnitude; if empty (default), use aliasing frequency
%   frac = half-length of the transition band (def=2 -> 1/2 octave). E.g.
%   if fc=623.89 Hz, smooth between 349.65 and 882.31 Hz. If frac==0, don't
%   smooth.
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
%   [1] Ben-Hur, Zamir, et al. "Efficient Representation and Sparse
%       Sampling of Head-Related Transfer Functions Using Phase-Correction
%       Based on Ear Alignment." IEEE/ACM Transactions on Audio, Speech,
%       and Language Processing 27.12 (2019): 2249-2262.
%   [2] Schörkhuber, C., Zaunschirm, M., & Höldrich, R. (2018). Binaural
%       Rendering of Ambisonic Signals via Magnitude Least Squares.
%       339–342.
%   [3] Zotter, Franz, and Matthias Frank. Ambisonics: A practical 3D audio
%       theory for recording, studio production, sound reinforcement, and
%       virtual reality. Springer Nature, 2019.
%
% AUTHOR: Isaac Engel (isaac.engel@imperial.ac.uk)
% February 2021

%% Some parameters
if ~exist('w','var')
    w = [];
end
if ~exist('r','var') || isempty(r)
    r = 0.085;
end
if ~exist('fc','var') || isempty(fc)
    c = 343; % speed of sound (m/s)
    fa = N*c/(2*pi*r); % aliasing frequency
    fc = max(3000,fa); % if fc not provided, use max of fa and 2kHz
end
if ~exist('frac','var') || isempty(frac)
    frac = 2;
end
if ~exist('earAz','var') || isempty(earAz)
    earAz = [pi/2, 3*pi/2];
end
if ~exist('earEl','var') || isempty(earEl)
    earEl = [pi/2, pi/2];
end
nfreqs = size(H,1);
f = linspace(0,fs/2,nfreqs).'; % frequency vector
c = 343; % speed of sound (m/s)
kr = 2*pi*f*r/c;
az = az(:).'; % force row vector
el = el(:).';

%% First, perform ear alignment
p = earAlign(kr,az,el,earAz,earEl);
H = H.*exp(-1i*p); % apply correction

%% Then, apply MagLS
[Hnm,fc] = toSH_MagLS(H,N,az,el,fs,w,fc,frac,r);

