function dtf = getDTF(h,fs)
% Get directional transfer function from an HRTF. Based on SOFAhrtf2dtf.
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021

flims = [50 20000]; % limits for CTF extraction
hlen = size(h,1);
nfreqs = hlen/2+1;
f = linspace(0,fs/2,nfreqs);
inds = f>=flims(1) & f<=flims(2);

H = ffth(h);
Hmag = abs(H);
CTFlogmag = mean(log(Hmag+eps),2); % mean of logmag rather than RMS
CTFlogmag_full = AKsingle2bothSidedSpectrum(CTFlogmag); % full spectrum
minphase = -imag( hilbert( CTFlogmag_full ) ); % minimum phase
CTF = exp(CTFlogmag).*exp(1i*minphase(1:nfreqs,:,:));
DTF = H;
DTF(inds,:,:) = H(inds,:,:)./CTF(inds,:,:);
ctf = iffth(CTF); % common transfer function
dtf = iffth(DTF); % directional transfer function
