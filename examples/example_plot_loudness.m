% Plot estimated loudness maps of two HRTFs, calculated with the 
% perceptually-weighted magnitude spectrum proposed in [1].
%
% EXTERNAL DEPENDENCIES:
%   SOFA API for Matlab (github.com/sofacoustics/API_MO)
%
% REFERENCES:
%   [1] Armstrong, Calum, et al. "A perceptual spectral difference model
%       for binaural signals." Audio Engineering Society Convention 145.
%       Audio Engineering Society, 2018.
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
nplotsH = 2; % number of plots in horizontal
nplotsV = 1; % number of plots in vertical
assert(n<=nplotsH*nplotsV,'Too few plots!')
clims = [];% colorbar limits, [] for automatic

%% Prepare figure
figsize = [16.4253 6.9003]; % change as appropriate
fig = figure('units','centimeters','pos',[2 2 figsize(1) figsize(2)],...
    'Renderer','painters','PaperSize',[figsize(1) figsize(2)],...
    'paperposition',[0 0 figsize(1) figsize(2)]);
gap = [.07 .008]; % gap between subplots in norm units (height width)
marg_h = [.13 .06]; % figure height margins in norm units (lower upper)
marg_w = [.06 .11]; % figure width margins in norm units (left right)
[ha, ~] = tight_subplot(nplotsV,nplotsH,gap,marg_h,marg_w);
colormap parula

%% Load HRTFs and plot loudness for all directions
for i=1:n
    SOFA_obj = SOFAload(filenames{i});
    [h,fs,az,el,r] = sofa2hrtf(SOFA_obj);
    H = ffth(h); % to frequency domain
    
    % Calculate loudness per direction
    [L,wERB] = perceptualSpectrum(H,fs);
    Lavgd = sum(mult2(L,wERB),1); % ERB-weighted avg loudness per direction

    axes(ha(i))
    plotSph(az,el,Lavgd(:,:,1)) % plot left ear only
    
    % Some plot settings
    title(labels{i}) % reference
    set(gca,'fontsize',7,'XDir','reverse')
    grid(gca,'off')
    if isempty(clims)
        clims = caxis;
    end
    caxis(clims)
    if i==n % colorbar only for the last plot
        c = colorbar; c.Label.String = 'Loudness (sones)';
        c.Position = [0.9 0.13 0.02 0.81]; % change as appropriate
    end
    if mod(i-1,nplotsH)==0 % show Y axis only for plots on the left
        yticks(-60:30:60)
        ylabel('Elevation (deg)')
    else
        yticklabels({})
    end
    if i>(nplotsH*(nplotsV-1)) % show X axis only for plots on last row
        xlabel('Azimuth (deg)')
    else
        xticklabels({})
    end
end
for i=(n+1):(nplotsH*nplotsV)
    axis(ha(i),'off') % hide unused axes
end
