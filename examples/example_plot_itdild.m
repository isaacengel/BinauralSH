% Compare ITDs and ILDs in polar plots for two HRTFs.
% ITD is calculated with MaxIACCe method, applying low-pass at 3 kHz [1]
% ILD is calculated according to [2]
%
% EXTERNAL DEPENDENCIES:
%   Auditory Modeling Toolbox (amtoolbox.sourceforge.net)
%
% REFERENCES:
%   [1] Andreopoulou, A., Katz, B.F.G., 2017. Identification of
%       perceptually relevant methods of inter-aural time difference
%       estimation. The Journal of the Acoustical Society of America 142,
%       588–598. https://doi.org/10.1121/1.4996457
%   [2] McKenzie, Thomas, Damian Murphy, and Gavin Kearney. "Interaural
%       level difference optimisation of first-order binaural Ambisonic
%       rendering." Audio Engineering Society Conference: 2019 AES
%       International Conference on Immersive and Interactive Audio. Audio
%       Engineering Society, 2019.
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021

%% Clear variables
clear

%% Parameters
basepath=which('binauralSH_start'); % base path
basepath=basepath(1:end-19); % Kill the function name from the path.

filenames = { % put here the paths to the HRTFs' SOFA files
    [basepath,'/hrtfs/FABIAN_HRIR_measured_HATO_0.sofa']
    [basepath,'/hrtfs/HRIR_L2702.sofa']
};
labels = { % put here the names for the legend
    'FABIAN'
    'KU100'
};
n = numel(filenames);
minEl = pi/180; % minimum elevation (in +/- radians) to be plotted

%% Prepare figure
figsize = [17.5472 7.3025];
fig = figure('units','centimeters','pos',[2 2 figsize(1) figsize(2)],...
    'Renderer','painters','PaperSize',[figsize(1) figsize(2)],...
    'paperposition',[0 0 figsize(1) figsize(2)]);
gap = [.04 .06]; % gap between subplots in norm units (height width)
marg_h = [0 .03]; % [.09 .1] figure height margins in norm units (lower upper)
marg_w = [.03 .23]; % figure width margins in norm units (left right)
[ha, ~] = tight_subplot(1,2,gap,marg_h,marg_w); % returns axis handles
pha = {[],[]}; % to put polar axis handles later
colors = parula(n+1);
ls = '-'; % line style
lw = 0.5; % line width
ms = 3; % marker size
mvec = '^vxsodph+'; % markers

%% Load HRTF and plot ITD and ILD for horizontal plane directions
for i=1:n
    SOFA_obj = SOFAload(filenames{i});
    [h,fs,az,el,r] = sofa2hrtf(SOFA_obj);
    % Find horizontal plane directions
    indHP = find(abs(el-pi/2)<minEl); % directions near the horz. plane
    azHP = az(indHP);
    % Calculate ITD and ILD for those directions
    itd = itdestimator(permute(h(:,indHP,:),[2,3,1]),'fs',fs,...
        'MaxIACCe','lp','upper_cutfreq', 3000,...
        'butterpoly', 10)*1e3; % requires AMT
    ild = getILD(h(:,indHP,:),fs);
    % Set how often we plot markers
    step_big = round(numel(azHP)/4);
    step_small = round(step_big/n);
    mi = ((i-1)*step_small+1):step_big:numel(azHP)-1; % marker indices
    % Plot ITD
    if isempty(pha{1})
        axes(ha(1))
    else
        axes(pha{1})
    end
    polarplot(azHP,abs(itd),'Color',colors(i,:),'LineWidth',lw,...
        'LineStyle',ls,'Marker',mvec(i),'MarkerSize',ms,'MarkerIndices',mi)
    pha{1} = gca;
    hold on
    % Plot ILD
    if isempty(pha{2})
        axes(ha(2))
    else
        axes(pha{2})
    end
    polarplot(azHP,abs(ild),'Color',colors(i,:),'LineWidth',lw,...
        'LineStyle',ls,'Marker',mvec(i),'MarkerSize',ms,'MarkerIndices',mi)
    pha{2} = gca;
    hold on
end

%% Change some plot settings
set(pha{1},'ThetaZeroLocation','top','ThetaLim',[0 360],'fontsize',7)
set(pha{1},'ThetaTick',[0:45:360],...
    'ThetaTickLabel',{'0','45','90','135','','-135','-90','-45'},...
    'RAxisLocation',-90)
title(pha{1},'ITD (ms)')
set(pha{2},'ThetaZeroLocation','top','ThetaLim',[0 360],'fontsize',7)
set(pha{2},'ThetaTick',[0:45:360],'ThetaTickLabel',{'0','45','90','135','','-135','-90','-45'},'RAxisLocation',-90)
title(pha{2},'ILD (dB)')
legend(labels,'position',[0.7986    0.3869    0.1991    0.2839]); % change size as appropriate