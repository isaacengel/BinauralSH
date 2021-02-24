function y = ffth(x,nfft,dim)
% Like fft(), but returns the spectrum only up to the Nyquist frequency.
%
% INPUT: same as fft()
%
% OUTPUT:
%   y = discrete Fourier transform of x, up to Nyquist frequency.
%
% EXTERNAL DEPENDENCIES:
%   AKtools (www.ak.tu-berlin.de/aktools)
%
% AUTHOR: Isaac Engel (isaac.engel@imperial.ac.uk)
% February 2021

if ~exist('nfft','var')
    nfft = [];
end

if ~exist('dim','var')
    dim = 1;
end

y = AKboth2singleSidedSpectrum(fft(x,nfft,dim));
    