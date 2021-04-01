function [hnm,fs,varOut] = toSH(h,N,varargin)
% Transform HRTF to SH domain at order N. This is a wrapper for the other
% toSH_* functions and applies some additional pre- and post- processing:
%   1. Zero-pad the HRIRs if indicated
%   2. Circular shift so the average HRIR onset falls at t=0
%   3. Obtain HRTF's SH coefficients using one of the available methods:
%       a. 'Trunc' = direct order truncation (default)
%       b. 'SpSub' = spatial subsampling (toSH_SpSub) [1]
%       c. 'MagLS' = magnitude least squares (toSH_MagLS) [2]
%       d. 'FDTA' = frequency-dependent time alignment (toSH_FDTA) [3]
%       e. 'TA' = phase correction by ear alignment (toSH_TA) [4]
%       f. 'SpSubTA' = frequency-dependent time alignment + spSub
%           (toSH_SpSubMod) [5]
%       g. 'BiMagLS' = phase correction by ear alignment + magLS
%   4. Apply tapering weights (e.g. Hann [6] or Max-rE [7])
%   5. Apply diffuse field EQ [8] or covariance constraint [9, sec.4.11.3]
%   6. Undo circular shift from step 2
%   7. Crop HRIR and apply fade in/out if indicated
%
% SIMPLE USAGE EXAMPLES:
%   hnm = toSH(h,15,'mode','MagLS','az',az,'el',el,'fs',48000);
%   hnm = toSH(sofa_h,15,'mode','MagLS');
%
% INPUT:
%   h = HRIRs in SOFA or matrix format (irlen x ndirs x 2 ears)
%   N = target SH order
%
% Mandatory inputs (as argument pairs) if h is provided in matrix format:
%   az = HRIR azimuth (ndirs x 1) in rad (0=front, pi/2=left)
%   el = HRIR elevation (ndirs x 1) in rad (0=top, pi/2=front)
%   fs = sampling frequency in Hz
%
% Optional inputs (provided as argument pairs):
%   mode = any of 'Trunc' (def), 'SpSub', 'MagLS', 'FDTA', 'TA', 
%       'SpSubMod', 'BiMagLS' (case insensitive; descriptions above)
%   w = quadrature weights (ndirs x 1); if empty (def), use pseudoinverse
%   EQ = one of the available equalisation options:
%       0 = no EQ (default)
%       1 = covariance constraint [9, sec.4.11.3]
%       2 = diffuse field EQ (HRF from [8])
%       3 = spherical head filters (SHF from [8]), adapted to tapering 
%           weights if appropriate [6]
%   fc = cutoff frequency for MagLS, frequency-dependent TA and dual-band 
%       tapering, in Hz (def=aliasing frequency)
%   nfftIn = if non-empty, zero-pad HRIRs to this length (samples)
%   nfftOut = if non-empty, crop the output to this length (samples)
%   r = head radius in m (def=0.085)
%   c = speed of sound in m/s (def=343)
%   earAz = left/right ear azimuth (1 x 2) in rad (def=[pi/2, 3*pi/2])
%   earEl = left/right ear elevation (1 x 2) in rad (def = [pi/2, pi/2])
%   Nmax = reference SH order for spSub and cov. constraint (def=35).
%   covConstPlots = whether to draw plots when processing covariance
%       constraint (def=false).
%   fadeIn = fade-in length in samples to apply to the output (def=0)
%   fadeOut = fade-out length in samples to apply to the output (def=0)
%   frac = half-length of the MagLS transition band (def=2 -> 1/2 octave).
%       E.g. if fc=623.89 Hz, smooth between 349.65 and 882.31 Hz. If
%       frac==0, don't smooth.
%   tapering = apply weights to high SH orders to reduce side lobes:
%       0 = off (default)
%       1 = original Hann weights from Hold et al (2019)
%       2 = improved (shorter) Hann window 
%       3 = max-rE weights from Zotter and Frank (2012), non-normalized
%       4 = same as 3, but normalize weights (McKenzie et al, 2018)
%   dualBand = whether to use dual-band tapering (def=false)
%   Hnm_ref = precomputed high-order Hnm used to calculate diffuse field
%       EQ. To speed things up if processing many HRTFs
%
% OUTPUTS:
%   hnm = HRIRs' SH coefficients (nfftOut x (N+1)^2 x 2 ears)
%   varOut = variable output, i.e. spherical head filter for 'Trunc', fc
%       for 'MagLS'/'BiMagLS', EQ filters if applied, etc.
%
% EXTERNAL DEPENDENCIES:
%   SOFA API for Matlab (github.com/sofacoustics/API_MO)
%   (only if the input HRTF is provided in SOFA format)
%
% REFERENCES:
%   [1] Bernschütz, Benjamin, et al. "Binaural reproduction of plane waves
%       with reduced modal order." Acta Acustica united with Acustica 100.5
%       (2014): 972-983.
%   [2] Schörkhuber, C., Zaunschirm, M., Höldrich, R., 2018. Binaural
%       Rendering of Ambisonic Signals via Magnitude Least Squares.
%       Presented at the DAGA 2018, Munich, Germany, pp. 339–342.
%   [3] Zaunschirm, M., Schörkhuber, C., Höldrich, R., 2018. Binaural
%       rendering of Ambisonic signals by head-related impulse response
%       time alignment and a diffuseness constraint. The Journal of the
%       Acoustical Society of America 143, 3616–3627.
%       https://doi.org/10.1121/1.5040489
%   [4] Ben-Hur, Zamir, et al. "Efficient Representation and Sparse
%       Sampling of Head-Related Transfer Functions Using Phase-Correction
%       Based on Ear Alignment." IEEE/ACM Transactions on Audio, Speech,
%       and Language Processing 27.12 (2019): 2249-2262.
%   [5] McKenzie, T., Murphy, D., Kearney, G., 2019. An evaluation of
%       pre-processing techniques for virtual loudspeaker binaural
%       ambisonic rendering, in: EAA Spatial Audio Signal Processing
%       Symposium. Paris, France, pp. 149–154.
%       https://doi.org/10.25836/sasp.2019.09
%   [6] Hold, C., Gamper, H., Pulkki, V., Raghuvanshi, N., Tashev, I.J.,
%       2019. Improving Binaural Ambisonics Decoding by Spherical Harmonics
%       Domain Tapering and Coloration Compensation, in: ICASSP 2019 - 2019
%       IEEE International Conference on Acoustics, Speech and Signal
%       Processing (ICASSP). Presented at the ICASSP 2019 - 2019 IEEE
%       International Conference on Acoustics, Speech and Signal Processing
%       (ICASSP), pp. 261–265. https://doi.org/10.1109/ICASSP.2019.8683751
%   [7] Zotter, Franz, and Matthias Frank. "All-round ambisonic panning and
%       decoding." Journal of the audio engineering society 60.10 (2012):
%       807-820.
%   [8] Ben-Hur, Zamir, et al. "Spectral equalization in binaural signals
%       represented by order-truncated spherical harmonics." The Journal of
%       the Acoustical Society of America 141.6 (2017): 4087-4096.
%   [9] Zotter, F., Frank, M., 2019. Ambisonics: A Practical 3D Audio
%       Theory for Recording, Studio Production, Sound Reinforcement, and
%       Virtual Reality, Springer Topics in Signal Processing. Springer
%       International Publishing, Cham.
%   [10] McKenzie, Thomas, Damian T. Murphy, and Gavin Kearney.
%       "Diffuse-field equalisation of binaural ambisonic rendering."
%       Applied Sciences 8.10 (2018): 1956.
%   [11] Archontis Politis, Microphone array processing for parametric
%       spatial audio techniques, 2016 Doctoral Dissertation, Department of
%       Signal Processing and Acoustics, Aalto University, Finland.
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021

% TODO: implement Tikhonov regularization for irregular grids
% According to Duraiswami 2004, it goes like this:
% reg_eps = 1e-6; % epsilon
% D = (1 + N*(N+1)) * eye((N+1)^2); % Tikhonov regularization matrix
% Y_inv = Y' * (Y * Y' + reg_eps*D)^(-1); % this instead of pseudoinverse

%% Parse inputs
p = inputParser;
addParameter(p,'mode','Trunc',@ischar) 
addParameter(p,'az',[])
addParameter(p,'el',[])
addParameter(p,'fs',[],@isscalar)
addParameter(p,'w',[])
addParameter(p,'EQ',0)
addParameter(p,'fc',[],@isscalar) 
addParameter(p,'nfftIn',[],@isscalar) 
addParameter(p,'nfftOut',[],@isscalar) 
addParameter(p,'r',0.0875,@isscalar)
addParameter(p,'c',343,@isscalar)
addParameter(p,'earAz',[pi/2, 3*pi/2])  
addParameter(p,'earEl',[pi/2, pi/2]) 
addParameter(p,'Nmax',35,@isscalar) 
addParameter(p,'covConstPlots',false,@isscalar) 
addParameter(p,'fadeIn',0,@isscalar);
addParameter(p,'fadeOut',0,@isscalar);
addParameter(p,'frac',2,@isscalar);
addParameter(p,'dualBand',false,@isscalar)
addParameter(p,'tapering',0,@isscalar)
addParameter(p,'Hnm_ref',[])

parse(p,varargin{:})
mode = p.Results.mode;
az = p.Results.az;
el = p.Results.el;
fs = p.Results.fs;
w = p.Results.w;
EQ = p.Results.EQ;
fc = p.Results.fc;
nfftIn = p.Results.nfftIn;
nfftOut = p.Results.nfftOut;
r = p.Results.r;
c = p.Results.c;
earAz = p.Results.earAz;
earEl = p.Results.earEl;
Nmax = p.Results.Nmax;
covConstPlots = p.Results.covConstPlots;
fadeIn = p.Results.fadeIn;
fadeOut = p.Results.fadeOut;
frac = p.Results.frac;
dualBand = p.Results.dualBand;
tapering = p.Results.tapering;
Hnm_ref = p.Results.Hnm_ref;

if isempty(fc) && ~strcmpi(mode,'BiMagLS')
    fc = N*c/(2*pi*r); % if fc not provided, use aliasing frequency
end

clear p
varOut = struct();

%% Process SOFA input
if isstruct(h) % if h is a struct, assume it is a SOFA object
    if ~isempty(fs) || ~isempty(az) || ~isempty(el)
        warning('Loading from SOFA overwrites parameters fs, az and el')
    end
    [h,fs,az,el] = sofa2hrtf(h); % overwrite az, el and fs
end

%% Preprocessing
% Zero-pad if indicated
if ~isempty(nfftIn)
    % Window
    win = hann(2*fadeIn-1);
    win = [ones(size(h,1)-fadeIn,1);win(fadeIn:end)];
    rms_before = rms(h);
    h = h.*win; % rename to h later
    rms_after = rms(h);
    nrg_loss = 100*(1-sum(rms_after(:))/sum(rms_before(:)));
    if nrg_loss > 1
        warning('%0.2f%% of the HRIR energy was lost when windowing. Consider disabling zero-padding or circ-shifting the HRIR to the right',nrg_loss)
    end
    % Zero pad
    h(end+1:nfftIn,:,:) = 0;
end

% Some parameters
irlen = size(h,1); 
nfft = 2^nextpow2(irlen);
nfreqs = nfft/2+1;
f = linspace(0,fs/2,nfreqs).'; % frequency vector
az = az(:); % force column vector
el = el(:);
w = w(:);

% Calculate delay so the SH signal is centered at t=0
% We do this by calculating the phase delay of the omni channel
% This is recommended for magLS calculations and doesn't hurt the other 
% methods

alignOption = 1; % 1=phase delay, 2=onset detection, 0=nothing
if alignOption == 1 % Option 1: phase delay 
    H = ffth(h,nfft,1); % to frequency domain
    Y0 = AKsh(0, [], az*180/pi, el*180/pi, 'real').';
    if ~isempty(w)
        Y0_inv = 4*pi*w.*Y0'; % if integrations weights are provided, use them
    else
        Y0_inv = pinv(Y0); % if not, the pseudoinverse will do just fine
    end
    H0 = mult3(H,Y0_inv); % 0th order SH signal
    pd = -unwrap(angle(H0(2:end,1,:)))./(2*pi*f(2:end)); % phase delay in s
    pd_mean = mean(pd,'all'); % avg delay
    H = H.*exp(1i*2*pi*f*pd_mean); % subtract delay
elseif alignOption == 2 % Option 2: onset detection
    onsL = AKonsetDetect(h(:,:,1)); % left ear onsets
    onsR = AKonsetDetect(h(:,:,2)); % right rear onsets
    ons = min([onsL;onsR]); % ipsilateral ear onsets
    %     nshift = -floor(min(ons))); % detect earliest onset
    nshift = -round(median(ons)); % median onset
    h = circshift(h,nshift);
    H = ffth(h,nfft,1);
else
    H = ffth(h,nfft,1);
end

%% Preprocess using the indicated method

if strcmpi(mode,'Trunc') % just order truncation
    % Order-truncated SH-HRTF (faster than order Nmax + truncation)
    Y = AKsh(N,[],az*180/pi,el*180/pi,'real').';
    if ~isempty(w)
        Y_inv = 4*pi*w.*Y'; % use integrations weights if provided
    else
        Y_inv = pinv(Y); % if not, pseudoinverse will do just fine
    end
    Hnm = mult3(H,Y_inv);

elseif strcmpi(mode,'SpSub') % subsampling aka virtual loudspeaker
    Hnm = toSH_SpSub(H,N,az,el,w,Nmax);

elseif strcmpi(mode,'FDTA') % frequency-dependent time-alignment
    if fc>=fs/2
        warning('fc (%0.2f Hz) >= fs/2 (%0.2f Hz). FDTA preprocessing has no effect...',fc,fs/2)
    end
    [Hnm,varOut.fc] = toSH_FDTA(H,N,az,el,fs,w,fc,r,earAz,earEl);

elseif strcmpi(mode,'MagLS') % magnitude least squares
    if fc>=fs/2
        warning('fc (%0.2f Hz) >= fs/2 (%0.2f Hz). MagLS preprocessing has no effect...',fc,fs/2)
    end
    [Hnm,varOut.fc] = toSH_MagLS(H,N,az,el,fs,w,fc,frac,r);

elseif strcmpi(mode,'TA') % time alignment (ear alignment)
    Hnm = toSH_TA(H,N,az,el,fs,w,r,earAz,earEl);

elseif strcmpi(mode,'BiMagLS') % ear alignment + MagLS
    if fc>=fs/2
        warning('fc (%0.2f Hz) >= fs/2 (%0.2f Hz). BiMagLS preprocessing has no effect...',fc,fs/2)
    end
    [Hnm,varOut.fc] = toSH_BiMagLS(H,N,az,el,fs,w,fc,frac,r,earAz,earEl);
        
elseif strcmpi(mode,'SpSubMod') % FDTA + SpSub
    if fc>=fs/2
        warning('fc (%0.2f Hz) >= fs/2 (%0.2f Hz). SpSubMod preprocessing is the same as SpSub...',fc,fs/2)
    end
    [Hnm,varOut.fc] = toSH_SpSubMod(H,N,az,el,fs,w,fc,r,earAz,earEl,Nmax);

else
    error('Unknown method %s. Use one of these: ''Trunc'', ''SpSub'', ''FDTA'', ''MagLS'', ''TA'', ''BiMagLS'', ''SpSubMod''',mode)
end

%% Tapering/max-rE weights

if tapering > 0
    
    if tapering == 4 && ~dualBand && EQ == 0
        warning('Normalized Max-rE weights (tapering=4) should be used with options dualBand=1 or EQ>0. Expect magnitude errors...')
    end
                    
    [~,wnm] = getTaperWin(N,tapering);

    if ~dualBand % apply tapering to whole signal [6]

        Hnm = Hnm.*wnm.';

    else % apply only to high frequencies [10]
        
        if fc>=fs/2
            warning('fc (%0.2f Hz) >= fs/2 (%0.2f Hz). Skipping dual band tapering...',fc,fs/2)
        else
            % Make crossover filters at fc (usually, f alias)
            % Filters are phase-matched and their sum ~= dirac delta
            lFilt = 100; % order of fir filters
            lpf = fir1(lFilt, fc/(fs/2), 'low').';
            hpf = fir1(lFilt, fc/(fs/2), 'high').';

            % Apply filters to the SH-HRTF
            nfreqs = size(H,1);
            Hnm_lf = Hnm.*ffth(lpf,2*(nfreqs-1)); % low-passed version
            Hnm_hf = Hnm.*ffth(hpf,2*(nfreqs-1)); % high-passed version
            p = linspace(0,pi*(lFilt/2+1),nfreqs).'; % filter phase delay
            Hnm_lf = Hnm_lf.*exp(1i*p); % undo filter delay
            Hnm_hf = Hnm_hf.*exp(1i*p); % undo filter delay

            % Apply weights to high frequency version
            Hnm_hf = Hnm_hf .* wnm.';

            % Output is the sum of the two
            Hnm = Hnm_lf + Hnm_hf;
        end
    end
    
end

%% Equalisation/covariance constraint

G1 = ones(nfreqs,2); % global EQ (left, right) 
G2 = zeros(nfreqs,2); % left/right mix (for covariance constraint)

if EQ > 0
    
    if EQ==1 || EQ==2 % these modes require a high-order reference
        if isempty(Hnm_ref)
            Y = AKsh(Nmax,[],az*180/pi,el*180/pi,'real').';
            if contains(mode,'ear') % time-align reference if required
                kr = 2*pi*f*r/c;
                p = earAlign(kr,az,el,earAz,earEl);
                Hear = H.*exp(-1i*p);
                Hnm_ref = mult3(abs(Hear),pinv(Y));
            else 
                Hnm_ref = mult3(abs(H),pinv(Y));
            end
        end
    end
    
    if EQ == 1 % covariance constraint [9, sec.4.11.3]
        
        [~,G1,G2] = applyCovConst(Hnm,Hnm_ref,covConstPlots,fs);
        
    elseif EQ == 2 % diffuse field EQ (HRF from [8])
        
        DF_high = sqrt( sum(abs(Hnm_ref).^2,2) )/(4*pi); % [8] eq.13
        DF_high = mean(DF_high,3); % average across ears    
        DF_N = sqrt( sum(abs(Hnm).^2,2) )/(4*pi);
        DF_N = mean(DF_N,3); % average across ears    
        % Option 1: direct inverse (could produce ringing artifacts)
%         G1 = h_high ./ h_N;
        % Option 2: regularized inverse (better)
        G1 = autoreg_minphase(DF_N,DF_high,fs);
        G1 = repmat(G1,1,2); % same for both ears
        
    elseif EQ == 3 % spherical head filters (SHF from [8])
        
        if ~strcmpi(mode,'Trunc')
            warning('SHF equalisation should only be used for non-preprocessed HRTFs; consider using EQ=3 (diffuse field EQ) instead')
        end
        if dualBand
            warning('SHF equalisation not optimised for dual band processing; consider using EQ=3 (diffuse field EQ) instead')
        end

        kr = 2*pi*f*r/c; 
        G1 = getSHF(N,kr,tapering);
        G1 = repmat(G1,1,2); % same for both ears
        
    end
    
    Hnm_EQ(:,:,1) = Hnm(:,:,1).*G1(:,1) + Hnm(:,:,2).*G2(:,1);
    Hnm_EQ(:,:,2) = Hnm(:,:,2).*G1(:,2) + Hnm(:,:,1).*G2(:,2);
    Hnm = Hnm_EQ;
    clear Hnm_EQ
 
end

varOut.G1 = G1;
varOut.G2 = G2;

%% Postprocessing
% Undo alignment at t=0
if alignOption == 1 % Option 1: phase delay
    Hnm = Hnm.*exp(-1i*2*pi*f*pd_mean);
    hnm = iffth(Hnm,[],1);
elseif alignOption == 2 % Option 2: onset detection + time shift
    hnm = iffth(Hnm,[],1);
    hnm = circshift(hnm,-nshift);
else
    hnm = iffth(Hnm,[],1);
end

% Crop/window output
rms_before = rms(hnm);
if nfftOut < nfft
    hnm(nfftOut+1:end,:,:) = []; % crop
end
fadeInVec = sin(linspace(0,pi/2,fadeIn).').^2; % fade in
fadeOutVec = cos(linspace(0,pi/2,fadeOut).').^2; % fade out
win = [fadeInVec;ones(size(hnm,1)-fadeIn-fadeOut,1);fadeOutVec]; % window
hnm = hnm.*win; % apply window
if nfftOut > nfft
    hnm(end+1:nfftOut,:,:) = 0; % zero-pad
end
rms_after = rms(hnm);
nrg_loss = 100*(1-sum(rms_after(:))/sum(rms_before(:)));
if nrg_loss > 1
    warning('%0.2f%% of the HRIR energy was lost when windowing. Consider changing fade in/out or output length.',nrg_loss)
end

