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
for ii = 1:size(CTFlogmag_full, 3)
  minphase(:,:,ii) = -imag( hilbert( CTFlogmag_full(:,:,ii) ) ); % minimum phase
end
CTF = exp(CTFlogmag).*exp(1i*minphase(1:nfreqs,:,:));
DTF = H;
DTF(inds,:,:) = H(inds,:,:)./repmat(CTF(inds,:,:), 1, size(H, 2), 1);
ctf = iffth(CTF); % common transfer function
dtf = iffth(DTF); % directional transfer function
