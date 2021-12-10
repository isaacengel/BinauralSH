function [EQmp, EQ] = autoreg_minphase(H,T,fs,maxAmp,frac,flims)
% Regularized inverse of an impulse response [REFERENCE PENDING].
%
% IMPORTANT: this function does not deal with phase and the output is
% always a minimum phase filter.
%
% SIMPLE USAGE EXAMPLE:
%   T = ones(nfreqs,1); % target is a flat response
%   EQ = autoreg_minphase(H,T,fs);
%
% INPUT:
%   H = original frequency response (nfreqs x nchannels)
%   T = target frequency response (nfreqs x 1)
%   fs = sampling rate in Hz
%   maxAmp = maximum amplification in dB (def=20)
%   frac = smoothing factor for notch inversion (def=4)
%   flims = inversion limits in Hz (def=[100 18000])
%
% OUTPUT:
%   EQmp = minimum phase EQ filter (nfreqs x nchannels)
%   EQ = linear phase EQ filter (nfrqes x nchannels)
%
% REFERENCES:
%   TODO
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021

%% Check inputs
if ~exist('maxAmp','var')
    maxAmp = 20; % default maximum amplification
end
if ~exist('frac','var')
    frac = 4; % default smoothing for sigma
end
if ~exist('flims','var')
    flims = [100,18000]; % default inversion limits
end

%% Calculate the inverse filter
X = H./T; % this is what we want to invert
Xmag = abs(X);
nfreqs = size(H,1);
f = linspace(0,fs/2,nfreqs).';

alpha = 1/(2*db2mag(maxAmp))^2;

Xhat = AKfractOctSmooth(Xmag,'amp',fs,frac,false);
sigma = zeros(nfreqs,1);
sigma(Xhat>=Xmag) = Xhat(Xhat>=Xmag)-Xmag(Xhat>=Xmag);

EQ = conj(X) ./ ( abs(X).^2 + alpha + sigma.^2 );
EQmag = abs(EQ);
EQph = angle(EQ);

%% Fade in/out to 1 outside the inversion limits
% Low-frequency boundary
fc1 = flims(1)*2^(-1/2); % smoothing starts 1/2 octaves before fc
fc2 = flims(1)*2^(1/2); % smoothing ends 1/2 octaves after fc
fc1_ind = find(f>fc1,1,'first'); % nearest available frequencies
fc2_ind = find(f<fc2,1,'last');
fc1 = f(fc1_ind);
fc2 = f(fc2_ind);
if ~isempty(fc1) && ~isempty(fc2) && (fc2-fc1)>0 
    fcvec_log = logspace(log10(fc1),log10(fc2),100); % f vector (log scale)
    w_log = linspace(0,1,100); % transition weights (log scale)
    fcvec = f(fc1_ind:fc2_ind); % f vector (lin scale)
    w = interp1(fcvec_log,w_log,fcvec,'spline'); % weights (lin scale)
    w = [zeros(fc1_ind-1,1);w(:);ones(nfreqs-fc2_ind,1)]; % for all freqs
    EQmag = (1-w)*1 + w.*EQmag; % weighted sum
end

% High-frequency boundary
fc1 = flims(2)*2^(-1/2); % smoothing starts 1/2 octaves before fc
fc2 = flims(2)*2^(1/2); % smoothing ends 1/2 octaves after fc
fc1_ind = find(f>fc1,1,'first'); % nearest available frequencies
fc2_ind = find(f<fc2,1,'last');
fc1 = f(fc1_ind);
fc2 = f(fc2_ind);
if ~isempty(fc1) && ~isempty(fc2) && (fc2-fc1)>0 
    fcvec_log = logspace(log10(fc1),log10(fc2),100); % f vector (log scale)
    w_log = linspace(1,0,100); % transition weights (log scale)
    fcvec = f(fc1_ind:fc2_ind); % f vector (lin scale)
    w = interp1(fcvec_log,w_log,fcvec,'spline'); % weights (lin scale)
    w = [ones(fc1_ind-1,1);w(:);zeros(nfreqs-fc2_ind,1)]; % for all freqs
    EQmag = (1-w)*1 + w.*EQmag; % weighted sum
end
EQ = EQmag.*exp(1i*EQph); % preserve original phase

%% Make minimum phase
EQmp = makeMinPhase(abs(EQ));

end
