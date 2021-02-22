function x = iffth(y,nfft,dim)
% Half-spectrum ifft

if ~exist('nfft','var')
    nfft = [];
end
if ~exist('dim','var')
    dim = 1;
end

x = ifft(AKsingle2bothSidedSpectrum(y),nfft,dim);