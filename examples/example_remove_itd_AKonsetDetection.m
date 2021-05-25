
basepath=which('binauralSH_start'); % base path
basepath=basepath(1:end-19); % Kill the function name from the path.

hrirname = [basepath,'/hrtfs/FABIAN_HRIR_measured_HATO_0.sofa'];

SOFA_obj = SOFAload(hrirname); % load HRTF in SOFA format
[h,fs,az,el] = sofa2hrtf(SOFA_obj); % get HRTF data

safety = 5; % number of safety samples

h_aligned = zeros(size(h));
for i=1:size(h,2)
    for j=1:size(h,3)
        hij = h(:,i,j);
        ons = AKonsetDetect(hij);
        ons = round(ons) - safety;
        hij(1:ons) = 0; % remove all data before onset
        hij = circshift(hij,-ons); % circular shift to the left
        h_aligned(:,i,j) = hij;
    end
end

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

