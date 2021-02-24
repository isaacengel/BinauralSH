function [Xson,w,calibrationGain] = perceptualSpectrum(X,fs,calibrationGain)
% Calculate perceptually-weighted magnitude spectrum according to [1].
% Return the ERB-weighted loudness in sons.
%
% INPUT:
%   X = frequency domain signal (half-spectrum) (nfreqs x nchannels x ears)
%   fs = sampling frequency in Hz
%   calibrationGain = optional input defining the gain to be applied to the
%       signal to calibrate it to db SPL (in dB; def=empty) 
%
% OUTPUT:
%   Xson = perceptually weighted spectrum in sones (nfreqs x nchannels x 
%       ears)
%   w = ERB weights for each frequency bin of Xson (nfreqs x 1)
%   calibrationGain = see above
%
% REFERENCES:
%   [1] Armstrong, Calum, et al. "A perceptual spectral difference model
%       for binaural signals." Audio Engineering Society Convention 145.
%       Audio Engineering Society, 2018.
%   [2] Glasberg, Brian R., and Brian CJ Moore. "Derivation of auditory
%       filter shapes from notched-noise data." Hearing research 47.1-2
%       (1990): 103-138.
%
% AUTHOR: Isaac Engel (isaac.engel@imperial.ac.uk)
% February 2021

if ~exist('calibrationGain','var')
    calibrationGain = [];
end

nfreqs = size(X,1);
f = linspace(0,fs/2,nfreqs).'; % frequency vector

Xdb = 20*log10(abs(X)); % signal to dB

% Get equal loudness contour (ISO 226:2003)
phonLvl = 60; % assuming 60 phons (60 dB SPL @ 1kHz ~= normal conversation)
[spl_sampled, f_sampled] = iso226(phonLvl);
f_sampled(end+1) = 22000; % add one last frequency to help interpolation
spl_sampled(end+1) = 100; % higher threshold, so interpolation goes up
spl = spline(f_sampled,spl_sampled,f); % interpolate
splToPhon = phonLvl - spl; % gain that converts dB SPL to Phon (0 at 1kHz)

% "Calibrate" the signal to dB SPL, given some calibration gain
if isempty(calibrationGain)
    % If a calibration value was not provided, normalize signal so the
    % average dB level at 1kHz is phonLvl (default = 60 dB SPL)
    [~,ind1khz] = min(abs(f-1000)); % find index for f=1kHz
    spl1khz = mean(Xdb(ind1khz,:,:),'all'); % average dB level at 1kHz
    calibrationGain = -spl1khz + phonLvl; % in dB
end
XdbSPL = Xdb + calibrationGain;

% Scale signals according to the contour (convert dB to phons)
Xphon = XdbSPL+splToPhon; 

% Convert phons to sones
Xson = 2.^((Xphon-40)/10);

% Apply ERB weighting (using formula from [2])
ERB = 21.4*log10(1+0.00437*f);
ERB(ERB==0) = ERB(2); % to prevent division by 0
w = 1./ERB; % weight = inverse of ERB
w = w./sum(w); % normalize so the sum of weights is 1
% ERB = ERB/sum(ERB); % normalize so the sum of weights is 1
% w = 1./ERB;
% Xavg = sum(Xson./ERB,1); % weighted average across frequencies

