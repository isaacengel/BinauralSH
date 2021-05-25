
basepath=which('binauralSH_start'); % base path
basepath=basepath(1:end-19); % Kill the function name from the path.

hrirname = [basepath,'/hrtfs/FABIAN_HRIR_measured_HATO_0.sofa'];

SOFA_obj = SOFAload(hrirname); % load HRTF in SOFA format
[h,fs,az,el] = sofa2hrtf(SOFA_obj); % get HRTF data

H = ffth(h); % to frequency domain
nfreqs = size(H,1); % number of frequency bins
f = linspace(0,fs/2,nfreqs).'; % frequency vector
r = 0.0875; % head radius (m)
c = 343; % speed of sound (m/s)
kr = 2*pi*f*r/c;

% Remove ITDs
p = earAlign(kr,az,el);
H = H.*exp(-1i*p);

h_aligned = iffth(H); % back to time domain

%% Plot before and after alignment for az=90, el=0
idx=SOFAfind(sofaobj,90,0); 
figure
subplot(1,2,1)
plot(squeeze(h(:,idx,:)))
AKp([h(:,idx,1),h(:,idx,2)],'t2d','fs',fs)
title('Original')
subplot(1,2,2)
AKp([h_aligned(:,idx,1),h_aligned(:,idx,2)],'t2d','fs',fs)
title('After removing ITD')
legend('Left HRIR', 'Right HRIR')
sgtitle('HRIRs for az=90, el=0')

