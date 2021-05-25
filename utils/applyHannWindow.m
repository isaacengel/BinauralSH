function y = applyHannWindow(x,fs,winlen_end,winlen_start)
% Apply Hann fade-in and/or fade-out to a time domain signal.
%
% SIMPLE USAGE EXAMPLE:
%   y = applyHannWindow(x,48000,0.01) % Hann window to last 10 ms of signal
%
% INPUT:
%   x = input signal in the time domain (nsamples x nchannels)
%   fs = sampling rate in Hz
%   winlen_end = fade-out length in seconds (def=0)
%   winlen_start = fade-in length in seconds (def=0)
%
% OUTPUT:
%   y = output signal (nsamples x nchannels)
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021

if ~exist('winlen_start','var')
    winlen_start = 0;
end
if ~exist('winlen_end','var')
    winlen_end = 0;
end

winlen_start = round(winlen_start*fs);
winlen_end = round(winlen_end*fs);

xlen = size(x,1);
nch = size(x,2);

if winlen_start > 0
    win1 = hann(2*winlen_start-1);
    win1 = [win1(1:winlen_start);ones(xlen-winlen_start,1)];
    win1 = repmat(win1,1,nch);
else
    win1 = ones(xlen,nch);
end

if winlen_end > 0
    win2 = hann(2*winlen_end-1);
    win2 = [ones(xlen-winlen_end,1);win2(winlen_end:end)];    
    win2 = repmat(win2,1,nch);
else
    win2 = ones(xlen,nch);
end

y = x .* win1 .* win2;