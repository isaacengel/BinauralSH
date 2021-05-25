function ild = getILD(h,fs)
% ILD claculation according to [1].
%
% INPUT:
%   h = HRIR (time x ndris x 2 ears)
%   fs = sampling frequency
%
% OUTPUT:
%   ild = ILD in dB
%
% REFERENCES:
%   [1] McKenzie, Thomas, Damian Murphy, and Gavin Kearney. "Interaural
%       level difference optimisation of first-order binaural Ambisonic
%       rendering." Audio Engineering Society Conference: 2019 AES
%       International Conference on Immersive and Interactive Audio. Audio
%       Engineering Society, 2019.
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021

hl = h(:,:,1);
hr = h(:,:,2);

% Linear phase high-pass FIR filters
fc = 1200; 
lFilt = 128; % order of fir filters
hpf = fir1(lFilt, fc/(fs/2), 'high').';

% Apply filter
nfft = max(4096,size(h,1));
Hl = ffth(hl,nfft);
Hr = ffth(hr,nfft);
Hpf = ffth(hpf,nfft);
%Hl_hp1 = Hl .* Hpf;
%Hr_hp1 = Hr .* Hpf;
Hl_hp = mult2(Hl,Hpf);
Hr_hp = mult2(Hr,Hpf);

% Get magnitude difference
md = 20*log10(abs(Hl_hp)./abs(Hr_hp));
    
% ERB space
nbands = 30;
lowFreq = 20;
highFreq = 20000;
cf = AMT_ERBSpace(lowFreq, highFreq, nbands);
cf = flipud(cf); % order from low to high
cf(end+1) = highFreq; % last frequency

% Loop through the ERB frequencies and get the average md for each one
f = linspace(0,fs/2,nfft/2+1);
ndirs = size(h,2);
ild_erb = zeros(nbands,ndirs);
for i=1:nbands
    ind1 = find(f>cf(i),1,'first');
    ind2 = find(f<=cf(i+1),1,'last');
    ild_erb(i,:) = mean(md(ind1:ind2,:),1);
end

% Average across ERB bands
ild = mean(ild_erb,1);
% 
% 
% fcoefs = MakeERBFilters(fs,nbands,lowFreq);
% 
% % Loop through all directions and get the ILD
% ndirs = size(h,2);
% ild = zeros(ndirs,1);
% for i=1:ndirs
%     % Apply filterbank and get magnitude difference
%     hl_hp_erb = ERBFilterBank(hl_hp(:,i), fcoefs);
%     hr_hp_erb = ERBFilterBank(hr_hp(:,i), fcoefs);
%     dif_erb = 20*log10(abs(hl_hp_erb)./abs(hr_hp_erb));
%     % Average across each band
%     dif = mean(dif_erb,
%     
% 
% %
% ild = 20*log10(abs(Hl_hp)./abs(Hr_hp));
%     % Apply filters to the subsampled HRTF
%     nfreqs = size(H,1);
%     Hs_lf = Hs.*ffth(lpf,2*(nfreqs-1)); % low-passed version
%     Hs_hf = Hs.*ffth(hpf,2*(nfreqs-1)); % high-passed version
%     p = linspace(0,pi*(lFilt/2+1),nfreqs).'; % phase delay of the filters
%     Hs_lf = Hs_lf.*exp(1i*p); % undo filter delay
%     Hs_hf = Hs_hf.*exp(1i*p); % undo filter delay
%     
% 
% onsL = AKonsetDetect(hl,10,-6,'rel',[3000,fs]);
% onsR = AKonsetDetect(hr,10,-6,'rel',[3000,fs]);
% 
% itd = 1e6*(onsL-onsR)/fs;
% 
