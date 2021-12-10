function [Hnm,fc] = toSH_MagLS(H,N,az,el,fs,w,fc,k,r,reg_eps)
% Transform HRTF to SH domain at order N with the MagLS method [1],
% following the simplified approach from [2] sec. 4.11.2 and 4.11.3. The
% main difference of this implementation vs [1] is that the phase is
% calculated iteratively in a simple way, instead of with a least squares
% approach.
%
% Important: the method works better if the HRIRs are aligned at t=0. This
% is taken care of in the wrapper function 'toSH'.
%
% SIMPLE USAGE EXAMPLE:
%   Hnm = toSH_MagLS(H,15,az,el,48000);
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
%   [1] Schörkhuber, C., Zaunschirm, M., Höldrich, R., 2018. Binaural
%       Rendering of Ambisonic Signals via Magnitude Least Squares.
%       Presented at the DAGA 2018, Munich, Germany, pp. 339–342.
%   [2] Zotter, F., Frank, M., 2019. Ambisonics: A Practical 3D Audio
%       Theory for Recording, Studio Production, Sound Reinforcement, and
%       Virtual Reality, Springer Topics in Signal Processing. Springer
%       International Publishing, Cham.
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021

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
Y = getRealSHmatrix(N,az,el);
Y_inv = getYinv(N,az,el,w,reg_eps);

%% Apply simplified MagLS from [2] sec.4.11.2

% For frequencies above a given cutoff, target the magnitude of the
% actual HRTF and the phase of the previous frequency bin. This avoids
% rapid changes in phase and the result is good enough for our
% purposes.

% Smooth transition parameters
fc1 = fc*2^(-k/2); % smoothing starts k/2 octaves before fc
fc2 = fc*2^(k/2); % smoothing ends k/2 octaves after fc
fc1_ind = find(f>fc1,1,'first'); % nearest available frequency
fc2_ind = find(f<fc2,1,'last'); % nearest available frequency
range1 = 1:fc1_ind-1; % first interval, below the cutoff (SFT)
range2 = fc1_ind:fc2_ind; % second interval (smooth transition)
range3 = fc2_ind+1:nfreqs; % third interval, above the cutoff (MagLS)

% Allocate space
Hnm = zeros(nfreqs,(N+1)^2,2); 

% Below the low cutoff, calculate SH coefficients via SFT
Hnm(range1,:,:) = mult3(H(range1,:,:),Y_inv);

% In the transition interval, calculate both ways and do weighted average
for fi=range2
    % First, calculate via SFT
    Hnm_sft = mult3(H(fi,:,:),Y_inv);
    % Second, calculate via MagLS
    Hprev = mult3(Hnm(fi-1,:,:),Y); % estimated HRTF for previous fi
    mag = abs(H(fi,:,:)); % magnitude of actual HRTF
    phase = angle(Hprev); % phase of previous fi
    Htarg = mag.*exp(1i*phase); % MagLS target HRTF
    Hnm_mls = mult3(Htarg,Y_inv); % MagLS estimation of the SH-HRTF
    % Finally, do a weighted average of both
    alpha = 0.5 + k/log(2) * log(f(fi)/fc); % 0 for fc1, 1 for fc2
    Hnm(fi,:,:) = alpha*Hnm_mls + (1-alpha)*Hnm_sft;
end

% Above the high cutoff, calculate SH coefficients via MagLS
for fi=range3
    Hprev = mult3(Hnm(fi-1,:,:),Y); % estimated HRTF for previous fi
    mag = abs(H(fi,:,:)); % magnitude of actual HRTF
    phase = angle(Hprev); % phase of previous fi
    Htarg = mag.*exp(1i*phase); % MagLS target HRTF
    Hnm(fi,:,:) = mult3(Htarg,Y_inv); % MagLS estimation of the SH-HRTF
end
