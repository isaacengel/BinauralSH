% Plot SH spectrum of an HRTF, i.e. the amount of energy per order of the 
% HRTF in the spherical harmonics domain, as in [1].
%
% EXTERNAL DEPENDENCIES:
%   SOFA API for Matlab (github.com/sofacoustics/API_MO)
%
% REFERENCES:
%   [1] Ben-Hur, Zamir, et al. "Efficient Representation and Sparse
%       Sampling of Head-Related Transfer Functions Using Phase-Correction
%       Based on Ear Alignment." IEEE/ACM Transactions on Audio, Speech,
%       and Language Processing 27.12 (2019): 2249-2262.
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021

%% Clear variables
clear

%% Parameters
N = 44; % maximum SH order
filenames = { % put here the paths to the HRTFs' SOFA files
    '../hrtfs/FABIAN_HRIR_measured_HATO_0.sofa'
    '../hrtfs/HRIR_L2702.sofa'
};
labels = { % put here the names for the legend
    'FABIAN'
    'KU100'
};
n = numel(filenames);
nplotsH = 2; % number of plots in horizontal
nplotsV = 1; % number of plots in vertical
assert(n<=nplotsH*nplotsV,'Too few plots!')

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

%% Load HRTF and plot ITD and ILD for horizontal plane directions
for i=1:n
    SOFA_obj = SOFAload(filenames{i});
    [hnm,fs] = toSH(SOFA_obj,N); % transform to SH domain with default settings
    Hnm = ffth(hnm); % to frequency domain (up to Nyquist frequency)
    Hnm_left = Hnm(:,:,1); % keep left ear only
    % Plot 
    axes(ha(i))
    plotSHenergy(Hnm_left,fs)
    title(labels{i})
    if i==n % colorbar only for the last plot
        c = colorbar; c.Label.String = 'Energy (dB)';
        c.Position = [0.9 0.13 0.02 0.81]; % change as appropriate
    end
    if mod(i-1,nplotsH)~=0 % show Y axis only for plots on the left
        yticklabels({})
        ylabel('')
    end
    if i<=(nplotsH*(nplotsV-1)) % show X axis only for plots on last row
        xticklabels({})
        xlabel('')
    end
    set(gca,'fontsize',7)
end
