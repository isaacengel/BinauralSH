function Hmp = makeMinPhase(Hmag)
% Takes a half-spectrum magnitude response and returns the minimum phase
% approximation of that signal, also half-spectrum. We assume the signal
% (Hmag) has an odd number of samples.
%
% INPUT:
%   Hmag = magnitude response up to Nyquist frequency (nfreqs x nchannels)
%
% OUTPUT:
%   Hmp = minimum phase reconstruction of Hmag
%
% REFERENCES:
%   TODO
%
% AUTHOR: Isaac Engel (isaac.engel@imperial.ac.uk)
% February 2021

Hlogmag = log(Hmag+eps); % to log-magnitude
Hlogmag_full = AKsingle2bothSidedSpectrum(Hlogmag); % full spectrum
phase = -imag( hilbert( Hlogmag_full ) ); % minimum phase
Hmp = exp(Hlogmag).*exp(1i*phase(1:end/2+1,:,:));