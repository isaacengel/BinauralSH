function plotSHenergy(Hnm,fs)
% Plot SH spectrum as in [1], Fig. 2.
%
% INPUT:
%   Hnm = SH coefficients up to Nyquist frequency (nfreqs x (N+1)^2 x 
%       nchannels). Several HRTFs can be concatenated like cat(3,Hnm1,Hnm2)
%   fs = sampling rate in Hz
%
% REFERENCES:
%   [1] Ben-Hur, Zamir, et al. "Efficient Representation and Sparse
%       Sampling of Head-Related Transfer Functions Using Phase-Correction
%       Based on Ear Alignment." IEEE/ACM Transactions on Audio, Speech,
%       and Language Processing 27.12 (2019): 2249-2262.
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021

%% Some variables
logscale = 0; % plot frequency in log scale
interpolateN = 1; % interpolate between SH orders to smooth the image
parulacolor = 0; % use parula colormap

%% Process inputs
nfreq = size(Hnm,1); % ,2); n freq bins
nSH = size(Hnm,2);
nch = size(Hnm,3);
N = sqrt(nSH)-1;
Nmax = max(N,44); % max N to plot
nfft = 2*(nfreq-1); % Hnm is half spectrum
f = (0:nfft/2)*fs/nfft;

% prepare axes for each channel
if nch==1
    ax = gca; 
else
    for ch=1:nch
        ax(ch) = subplot(1,nch,ch); 
    end
end

%% Get SH spectrum
% Normalize so the spectra plot peaks at 0dB
nrg = sum(abs(Hnm).^2,2); % energy per channel
maxnrg = max(nrg(:)); % peak value
Hnm = Hnm./sqrt(maxnrg);

e = zeros(nfreq,N+1,nch); % energy per order
ecum = zeros(nfreq,N+1,nch); % energy accumulated until each order
for n=0:N
    if n==0
        ch = 1;
    else
        ch = (n^2+1):((n+1)^2);
    end
    e(:,n+1,:) = sum(abs(Hnm(:,ch,:)).^2,2);
    ecum(:,n+1,:) = sum(e,2);
end

edb = 10*log10(e);
ecumdb = 10*log10(ecum);

% warp to log scale if indicated
if logscale
    fvec = logspace(log10(f(2)),log10(f(end)),nfreq); % f vector (log scale)
    edb = interp1(f,edb,fvec,'spline');
    ecumdb = interp1(f,ecumdb,fvec,'spline');
else
    fvec = f;
end

% interpolate if indicated
if interpolateN
    npoints = 1000;
    shvec = linspace(0,Nmax,npoints);
    if N<Nmax
        ind_int = shvec <= N;
    else
        ind_int = 1:npoints;
    end
    edb_int = nan(nfreq,npoints,nch);
    for ch=1:nch
        edb_int(:,ind_int,ch) = interp1(0:N,edb(:,:,ch).',shvec(ind_int),'spline').';
        edb_int(:,~ind_int,ch) = -300; % very low value
    end
    edb = edb_int;
else
    shvec = 0:N;
end

%% Plot
for ch = 1:nch
    axes(ax(ch))

    h = imagesc(fvec,shvec,edb(:,:,ch).');

    % colorbar
    if parulacolor
        colormap(parula)
    else
        colormap(flipud(gray))
    end

    ylabel('SH order')

    xlim([0 20000]) % 100
    if logscale
        set(gca,'XScale','log')
        xticks([100,1000,10000,20000])
        xticklabels({'100','1k','10k','20k'})    
        xlabel('f (Hz)')
    else    
        xticks([1000,10000,20000])
        xticklabels({'1','10','20'})
        set(gca,'XScale','lin')
        xlabel('f (kHz)')
    end
    set(gca,'YDir','normal')
    set(gca,'clim',[-40 0]) % [-22.5 7.5])
    % fprintf('maxdB = %0.2f, mindB = %0.2f\n',max(edb(:)),min(edb(:)))
    view(2)
    ylim([0 Nmax])

    % Find b90 and b99
    b90 = 0;
    b99 = 0;
    %b90 = sum(~(ecum(:,:,ch) > 0.9*ecum(:,end,ch)),2);
    %b99 = sum(~(ecum(:,:,ch) > 0.99*ecum(:,end,ch)),2);
    for ii = 1:size(ecum, 2) % for compatibility with Matlab 2016
        b90 = b90 + ~(ecum(:,ii,ch) > 0.9*ecum(:,end,ch));
        b99 = b99 + ~(ecum(:,ii,ch) > 0.99*ecum(:,end,ch));
    end
    % Plot them
    hold on
    plot(fvec,b90,':','LineWidth',1,'Color','r')
    plot(fvec,b99,':','LineWidth',1,'Color',[0 0.75 0])
    legend({'b90','b99'},'location','nw')
    hold off
end
