function Y = smoothERB(X,fs,factor,n)
% Apply ERB smoothing to a magnitude spectrum following [1,2].
%
% SIMPLE USAGE EXAMPLE:
%   H = ffth(h); % spectrum up to Nyquist
%   Hsm = smoothERB(H,fs); % smoothed signal
%
% INPUT:
%   X = spectrum up to Nyquist frequency (nfreqs x nch)
%   fs = sampling rate in Hz
%   factor = size of the smoothing window in ERB (def=1)
%   n = filter order (def=4)
%
% OUTPUT:
%   Y = smoothed signal (nfreqs x nch)
%
% REFERENCES:
%   [1] Kohlrausch, Armin, and Jeroen Breebaart. "Perceptual (ir) relevance
%       of HRTF magnitude and phase spectra." Audio Engineering Society
%       Convention 110. Audio Engineering Society, 2001.
%   [2] Hassager, Henrik Gert, Fredrik Gran, and Torsten Dau. "The role of
%       spectral detail in the binaural transfer function on perceived
%       externalization in a reverberant environment." The Journal of the
%       Acoustical Society of America 139.5 (2016): 2992-3000.
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% May 2021

%% Some parameters
if ~exist('n','var')
    n = 4;
end
if ~exist('factor','var')
    factor = 1;
end
nfreqs = size(X,1); % number of frequency bins
f = linspace(0,fs/2,nfreqs).'; % vector with frequencies
Xmag2 = abs(X).^2; % |X(f)|^2
Y = zeros(size(X));

%% Smooth according to [1]
b = 24.7*(0.00437*f+1) ./ (2*sqrt(2^(1/n)-1)); % [1] Eq. 4 for all fc=f
b = b*factor; % scale ERB factor as in [2]
for ind=1:nfreqs
    fc = f(ind);
    H = ( 1+1i*(f-fc)./b(ind) ).^(-n); % Eq. 2
    Hmag2 = abs(H).^2; % |H(f,fc)|^2
    Y(ind,:,:) = sqrt( sum(Xmag2.*Hmag2,1) ./ sum(Hmag2,1) );
end