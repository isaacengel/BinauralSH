function G = getSHF(N,kr,option,Nhigh,doplot)
% Return spherical head filters (SHF) for order N, following [1].
% Optionally, adapt the filters to tapering weights, as in [2].
%
% NOTE: we get numerical errors for Nhigh>92
%
% INPUT:
%   N = target SH order
%   kr = kr vector (nfreqs x 1)
%   option = if > 0, adapt SHF to tapering weights (see getTaperWin)
%   Nhigh = reference SH order (def = 90). 
%   doplot = plot SHF (def = 0)
%
% OUTPUT:
%   G = spherical head filters (nfreqs x 1)
%
% EXTERNAL DEPENDENCIES:
%   AKtools (www.ak.tu-berlin.de/aktools)
%
% REFERENCES:
%   [1] Ben-Hur, Zamir, et al. "Spectral equalization in binaural signals
%       represented by order-truncated spherical harmonics." The Journal of
%       the Acoustical Society of America 141.6 (2017): 4087-4096.
%   [2] Hold, C., Gamper, H., Pulkki, V., Raghuvanshi, N., Tashev, I.J.,
%       2019. Improving Binaural Ambisonics Decoding by Spherical Harmonics
%       Domain Tapering and Coloration Compensation, in: ICASSP 2019 - 2019
%       IEEE International Conference on Acoustics, Speech and Signal
%       Processing (ICASSP). Presented at the ICASSP 2019 - 2019 IEEE
%       International Conference on Acoustics, Speech and Signal Processing
%       (ICASSP), pp. 261–265. https://doi.org/10.1109/ICASSP.2019.8683751
%   [3] Bernschütz, B., 2016. Microphone arrays and sound field
%       decomposition for dynamic binaural recording (Doctoral Thesis).
%       Technische Universität Berlin. Technische Universität Berlin,
%       Berlin. https://doi.org/10.14279/depositonce-5082
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021

%% Fill missing parameters
if ~exist('option','var')
    option = 0;
end
if ~exist('Nhigh','var')
    Nhigh = 90; % >=93 gives numerical errors
end
if ~exist('doplot','var')
    doplot = false;
end

%% Process 
nfreqs = numel(kr);
bn = zeros(Nhigh+1,nfreqs);

% Replace zeros in kr with the next value
for i=1:length(kr)
    if kr(i) <= 0
        for j=i+1:length(kr)
            if kr(j) >= 0
                kr(i) = kr(j);
                break
            end
        end
    end
end

% Pn for all orders
for n=0:Nhigh
    % bn filter calculation (see [3], eq.3.20)
    bessel1 = AKshRadial(kr,'bessel',1,n,false);
    dbessel1 = AKshRadial(kr,'bessel',1,n,true);
    hankel2 = AKshRadial(kr,'hankel',2,n,false);
    dhankel2 = AKshRadial(kr,'hankel',2,n,true);
    bn(n+1,:) = 4*pi*(1i)^n*(bessel1-(dbessel1./dhankel2).*hankel2);
end

nvec_N = diag(2*(0:N)+1);
if option > 0
    wN = diag(getTaperWin(N,option));
else
    wN = 1;
end
% Adapted weighs ([2], eq. 10). NOTE: wN is squared here, but not in the
% original paper (their mistake?)
pn_N = sqrt( sum( wN.^2*nvec_N*abs(bn(1:N+1,:)).^2 ,1) )/(4*pi); 

nvec_ref = diag(2*(0:Nhigh)+1);
pn_ref = sqrt( sum( nvec_ref*abs(bn(1:Nhigh+1,:)).^2 ,1) )/(4*pi);

% Resulting EQ filter
G = pn_ref./pn_N; % keep only the filter for order N
G = G.'; % transpose to column vector

if doplot
    figure, semilogx(kr,db(abs(G))), xlim([0.03 38]), ylim([0 30])
end

% Plots ([1] fig. 2-3)
% figure, semilogx(kr,db(abs(pn([4 7 13 31],:).')))
% xlim([0.03 38]), ylim([-20 5])
