function y = ffth(x,nfft,dim)
% Half-spectrum fft

if ~exist('nfft','var')
    nfft = [];
end
if ~exist('dim','var')
    dim = 1;
end

y = AKboth2singleSidedSpectrum(fft(x,nfft,dim));
    