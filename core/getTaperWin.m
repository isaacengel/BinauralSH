function [wN,wnm] = getTaperWin(N,option)
% Return tapering weights for order N.
%
% INPUT:
%   N = maximum SH order
%   option = choose type of tapering weights:
%       1 = Hann window according to [1]
%       2 = improved Hann window, applied only to 3 highest SH orders, to
%           prevent excessive attenuation
%       3 = Max-rE weights [2]
%       4 = normalized Max-rE weights [3]
%   tapering = whether to use the tapering version of SHFs (def=0),
%       1=Hanning window, 2=maxRE weights
%   Nhigh = reference SH order (def = 90). Note that we get numerical
%       errors for Nhigh>92
%   doplot = plot SHF
%
% OUTPUT:
%   G = spherical head filters ( nfreqs x 1 )
%
% REFERENCES:
%   [1] Hold, C., Gamper, H., Pulkki, V., Raghuvanshi, N., Tashev, I.J.,
%       2019. Improving Binaural Ambisonics Decoding by Spherical Harmonics
%       Domain Tapering and Coloration Compensation, in: ICASSP 2019 - 2019
%       IEEE International Conference on Acoustics, Speech and Signal
%       Processing (ICASSP). Presented at the ICASSP 2019 - 2019 IEEE
%       International Conference on Acoustics, Speech and Signal Processing
%       (ICASSP), pp. 261–265. https://doi.org/10.1109/ICASSP.2019.8683751
%   [2] Zotter, Franz, and Matthias Frank. "All-round ambisonic panning and
%       decoding." Journal of the audio engineering society 60.10 (2012):
%       807-820.
%   [3] McKenzie, Thomas, Damian T. Murphy, and Gavin Kearney.
%       "Diffuse-field equalisation of binaural ambisonic rendering."
%       Applied Sciences 8.10 (2018): 1956.
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021

%% Check input
if ~exist('option','var')
    option = 1;
end

%% Get tapering weights
if option == 1
    % Long Hann window [1]
    halflen_hann = floor((N + 1)/2) + 1;
    w_full = hann(2 * halflen_hann); % full hann window
    wN = ones(N + 1,1); % half hann window
    wN(end-halflen_hann+2:end) = w_full(end-halflen_hann+1:end-1);
    % vector with one index per nm
    nm = AKgetNM(N); % get nm indices
    wnm = zeros((N+1)^2,1);
    for n=0:N
        ind = nm == n;
        wnm(ind) = wN(n+1);
    end
elseif option == 2
    % Improved, shorter Hann window
    halflen_hann = min(4,floor((N + 1)/2) + 1); % top out at 4
    w_full = hann(2 * halflen_hann); % full hann window
    wN = ones(N + 1,1); % half hann window
    wN(end-halflen_hann+2:end) = w_full(end-halflen_hann+1:end-1);
    % vector with one index per nm
    nm = AKgetNM(N); % get nm indices
    wnm = zeros((N+1)^2,1);
    for n=0:N
        ind = nm == n;
        wnm(ind) = wN(n+1);
    end
elseif option == 3
    % Max-rE weights [2]
    wnm = getMaxREweights(N);
    wN = unique(wnm,'stable');
elseif option == 4
    % Normalized Max-rE weights [3]
    wnm = getMaxREweights(N);
    wnm = wnm./rms(wnm);
    wN = unique(wnm,'stable');
else
    error('wrong tapering option %d',option)
end
