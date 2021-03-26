function [Hnm,fc] = toSH_MagLS(H,N,az,el,fs,w,fc,frac,r)
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
%   frac = half-length of the transition band (def=2 -> 1/2 octave). E.g.
%       if fc=623.89 Hz, smooth between 349.65 and 882.31 Hz. If frac==0, 
%       don't smooth.
%   r = head radius in m (def=0.085)
%   
% OUTPUT:
%   Hnm = HRTF's SH coefficients (nfreqs x (N+1)^2 x 2 ears)
%   fc = see above
%
% EXTERNAL DEPENDENCIES:
%   AKtools (www.ak.tu-berlin.de/aktools)
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
    r = 0.085; % default head radius
end
if ~exist('fc','var') || isempty(fc)
    c = 343; % speed of sound (m/s)
    fc = N*c/(2*pi*r); % if fc not provided, use aliasing frequency
end
if ~exist('frac','var') || isempty(frac)
    frac = 2; % default smoothing = half octave
end
nfreqs = size(H,1);
f = linspace(0,fs/2,nfreqs).'; % frequency vector

%% Get the order-truncated SH-HRTF
Y = AKsh(N, [], az*180/pi, el*180/pi, 'real').'; 
if ~exist('w','var') || ~isempty(w)
    Y_inv = 4*pi*w.*Y'; % if integrations weights are provided, use them
else
    Y_inv = pinv(Y); % if not, the pseudoinverse will do just fine
end
Hnm = mult3(H,Y_inv); % SH coeffs up to N

%% Apply simplified MagLS from [2] sec.4.11.2

% For frequencies above a given cutoff, target the magnitude of the
% actual HRTF and the phase of the previous frequency bin. This avoids
% rapid changes in phase and the result is good enough for our
% purposes.

% Smooth transition parameters
if frac ~= 0
    % Option 1: fc is at the middle of the transition band
    fc1 = fc*2^(-1/frac); % smoothing starts 1/frac octaves before fc
    fc2 = fc*2^(1/frac); % smoothing ends 1/frac octaves after fc
    % Option 2: fc is at the end of the transition band
%     fc1 = fc*2^(-1/(2*frac));
%     fc2 = fc;
    % Look for nearest available frequencies
    fc1_ind = find(f>fc1,1,'first'); 
    fc2_ind = find(f<fc2,1,'last');
    fc1 = f(fc1_ind);
    fc2 = f(fc2_ind);
    % Define frequency and weight vector
    fcvec_log = logspace(log10(fc1),log10(fc2),100); % f vector (log scale)
    sw_log = linspace(0,1,100); % transition weights (log scale)
    fcvec = f(fc1_ind:fc2_ind); % f vector (lin scale)
    sw = interp1(fcvec_log,sw_log,fcvec,'spline'); % weights (lin scale)
    sw = [zeros(fc1_ind-1,1);sw(:);ones(nfreqs-fc2_ind,1)]; % for all freqs
else
    sw = ones(nfreqs,1);
    fc1_ind = find(f>fc,1,'first');
end

for fi=fc1_ind:nfreqs % for frequencies above the cutoff
    prevH = mult3(Hnm(fi-1,:,:),Y); % previous fi estimated HRTF
    mag = abs(H(fi,:,:)); % magnitude of actual HRTF
    phase = angle(prevH); % phase of previous fi
    currH = mag.*exp(1i*phase); % current fi estimated HRTF
    newHnm = mult3(currH,Y_inv); % magLS estimation of the SH-HRTF
    Hnm(fi,:,:) = sw(fi)*newHnm + (1-sw(fi))*Hnm(fi,:,:); % smoothing
end

