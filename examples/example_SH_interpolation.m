% Simple SH interpolation example
%
% EXTERNAL DEPENDENCIES:
%   SOFA API for Matlab (github.com/sofacoustics/API_MO)
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% March 2021

%% Clear variables
clear

%% Parameters
N = 44; % SH order (~40 is high enough for very accurate interpolation)
filename = '../hrtfs/FABIAN_HRIR_measured_HATO_0.sofa';

%% Load HRTF
SOFA_obj = SOFAload(filename); % load HRTF in SOFA format
[h,fs,az,el] = sofa2hrtf(SOFA_obj); % get HRTF data
ndirs = size(h,2); % number of directions

%% Remove a few directions to interpolate them later
alldirs = true(ndirs,1);
alldirs(randi(ndirs,100,1)) = false; % remove 100 random directions
targetdirs = ~alldirs; % this are the directions that we will interpolate

% Test directions
h_target = h(:,targetdirs,:);
az_target = az(targetdirs);
el_target = el(targetdirs);

% Remaining directions used to interpolate
h_toInt = h(:,alldirs,:);
az_toInt = az(alldirs);
el_toInt = el(alldirs);

%% Put HRTF in SH domain
hnm = toSH(h_toInt,N,'az',az_toInt,'el',el_toInt,'fs',fs);

%% Interpolate test directions
h_int = fromSH(hnm,fs,az_target,el_target);

%% Do some plots for one direction
figure('pos',[25.8000 53 930.4000 282.4000])
subplot(1,3,1)
AKp([h_target(:,1,1),h_int(:,1,1)],'t2d','fs',fs)
legend({'Target','Interpolated'},'location','northeast')
title('Time domain')
subplot(1,3,2)
AKp([h_target(:,1,1),h_int(:,1,1)],'m2d','fs',fs)
title('Magnitude')
subplot(1,3,3)
AKp([h_target(:,1,1),h_int(:,1,1)],'pu2d','fs',fs)
title('Unwrapped phase')
sgtitle(sprintf('Interpolated vs target (HRTF at az=%0.2f, el=%0.2f deg)',rad2deg(az_target(1)),rad2deg(el_target(1))))
