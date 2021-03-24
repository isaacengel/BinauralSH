% SH interpolation example. This script performs HRTF interpolation using
% two different methods:
%   1. Trunc (SH order truncation): baseline SH interpolation without
%       additional preprocessing.
%   2. TA (time-alignment): HRTFs are time-aligned before interpolation,
%       which improves the performance at low orders.
%
% For each method, we evaluate spatial orders ranging from 1 (lowest usable
% order) to 40 (very high order).
%
% The following plots are shown:
%   1. Interpolation errors (magnitude and phase) vs the target HRTF.
%   2. Interaural time differences (ITD) and level differences (ILD) and
%       their error vs the target HRTF.
%   3. Loudness maps and perceptual spectral difference (PSD).
%   4. Performance per order according to various binaural models:
%       - Reijniers2014: localisation performance (lateral and polar)
%       - Baumgartner2020: externalisation rates
%       - Jelfs2011: speech reception in presence of maskers
%
% WARNING: this script can take several minutes to run and will generate a
% few large files (~400MB in total)
%
% EXTERNAL DEPENDENCIES:
%   SOFA API for Matlab (github.com/sofacoustics/API_MO)
%   Auditory Modeling Toolbox (amtoolbox.sourceforge.net)
%   AKtools (www.ak.tu-berlin.de/aktools)
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% March 2021

%% Clear variables
clear

%% Parameters
N_vec = [1,5,10:10:40]; % spatial orders to test
Nref = 44; % reference (very high) spatial order
r = 0.0875; % head radius = 8.75 cm
c = 343; % speed of sound
itdJND = 20; % ITD JND for plots, according to Klockgether 2016
ildJND = 0.6; % ILD JND for plots, according to Klockgether 2016

% HRTF file
hrirname = '../hrtfs/FABIAN_HRIR_measured_HATO_0.sofa'; % change as needed

% Working directory
workdir = '../processed_data'; % change as needed
[~,hrtfname] = fileparts(hrirname);
workdir = [workdir,'/',hrtfname];
if ~isfolder(workdir)
    fprintf('Making folder %s...\n',workdir)
    mkdir(workdir);
end
if ~isfolder([workdir,'/hnm'])
    fprintf('Making folder %s...\n',[workdir,'/hnm'])
    mkdir([workdir,'/hnm'])
end
if ~isfolder([workdir,'/results'])
    fprintf('Making folder %s...\n',[workdir,'/results'])
    mkdir([workdir,'/results'])
end

%% Load HRTF
SOFA_obj = SOFAload(hrirname); % load HRTF in SOFA format
[h,fs,az,el] = sofa2hrtf(SOFA_obj); % get HRTF data

% Zero-pad the HRIRs to increase frequency resolution
taps = 2048;
h = [h;zeros(taps-size(h,1),size(h,2),size(h,3))];

H = ffth(h); % to frequency domain
nfreqs = size(H,1);
f = linspace(0,fs/2,nfreqs).';

%% Get high-order Hnm
filename = sprintf('%s/hnm/ref.mat',workdir);
if isfile(filename) % load if it already exists
    fprintf('Found %s. Loading...\n',filename)
    load(filename,'hnm_ref','Hnm_TA_ref')
else
    fprintf('Making reference hnm (also time-aligned Hnm, for Fig. 1)...\n');
    % Standard
    hnm_ref = toSH(h,Nref,'az',az,'el',el,'fs',fs,'mode','trunc');
    % Time-aligned
    hnm_TA_ref = toSH(h,Nref,'az',az,'el',el,'fs',fs,'mode','TA');
    Hnm_TA_ref = ffth(hnm_TA_ref);
    % Save
    fprintf('Saving %s...\n',filename)
    save(filename,'hnm_ref','Hnm_TA_ref','fs')
end
Hnm_ref = ffth(hnm_ref);

%% Defining test conditions
ncond = 2; 
% 1. "Trunc": baseline method, with no preprocessing such as time-alignment
test_conditions{1}.name = 'Trunc';
test_conditions{1}.preproc = 'trunc';
% 2. "TA": time-align the HRTF prior to getting the SH coefficients
test_conditions{2}.name = 'TA';
test_conditions{2}.preproc = 'TA';

%% Define a few direction subsets
% 1. Nearest neighbours to 110-point Lebedev grid
gridTestReq = sofia_lebedev(110,0);
[~,indLeb] = getGridSubset([az,el],gridTestReq,0);
indLeb = find(indLeb); % get indices
azLeb = az(indLeb);
elLeb = el(indLeb);
wLeb = gridTestReq(:,3);
% 2. Median plane (lat = 0 +/- 1 deg)
[lat,~]=sph2hor(az*180/pi,90-el*180/pi); % lat/pol coordinates
indMP = find(abs(lat)<1); 
azMP = az(indMP);
elMP = el(indMP);
% 3. Horizontal plane (el = 90 +/- 1 deg)
indHP = find(abs(el-pi/2)<(pi/180));
azHP = az(indHP);
% elHP = el(indHP);
indFront = find(abs(el-pi/2)<(pi/180) & abs(az)<(pi/180)); % front

%% Generate and save hnm for each test condition
fprintf('Generating hnms...\n')
for N=N_vec % iterate through spatial orders
    fprintf('Processing order %d...\n',N);
    for i=1:ncond
        name = test_conditions{i}.name;
        preproc = test_conditions{i}.preproc;
        filename = sprintf('%s/hnm/ord%0.2d_%s.mat',workdir,N,name);
        if isfile(filename)
            fprintf('\tFound %s. Skipping...\n',filename)
        else
            fprintf('\tGenerating hnm...\n');
            [hnm,varOut] = toSH(h,N,'az',az,'el',el,'fs',fs,'mode',preproc);
            fprintf('\tSaving %s...\n',filename)
            save(filename,'hnm','varOut')
        end  
    end
end

%% Interpolate HRTFs and generate results
fprintf('Generating results...\n')

%% First, run analysis for the reference (original HRTF)
filename = sprintf('%s/results/ref.mat',workdir);
if isfile(filename)
    fprintf('\tFound %s. Loading...\n',filename)
    load(filename,'results')
else
    fprintf('\tGenerating results for reference HRTF...\n')

    %% Numerical analysis
    % Magnitude and phase delay 110-point Lebedev grid
    results.mag = 20*log10(abs(H(:,indLeb,:)));
    results.pd = (-unwrap(angle(H(:,indLeb,1)))+unwrap(angle(H(:,indLeb,2))))./(2*pi*f)*1e6;
    % Loudness across directions
    [L,wERB,calibrationGain] = perceptualSpectrum(H,fs);
    Lavgd = sum(L.*wERB,1); % ERB-weighted avg loudness per direction
    results.L = Lavgd;
    results.L_perDirection = L; % save to calculate PSD
    results.Lavg = sum(Lavgd(:,indLeb,:).*wLeb.',2); % weighted avg on the Lebedev grid
    results.calibrationGain = calibrationGain; % save to calibrate loudness
    % Interaural differences for horizontal plane
    results.itd = itdestimator(permute(h(:,indHP,:),[2,3,1]),'fs',fs,...
        'MaxIACCe','lp','upper_cutfreq', 3000,...
        'butterpoly', 10)*1e6;
    results.ild = getILD(h(:,indHP,:),fs);

    %% Reijniers 2014
    fprintf('\t\tRunning reijniers2014...\n'); tic
    % Make DTF and SOFA object for Lebedev grid directions only
    dtf = getDTF(h(:,indLeb,:),fs);
    SOFA_obj = makesofa(dtf,fs,azLeb,elLeb);
    % Preprocessing source information (demo_reijniers2014)
    [template_loc, target] = reijniers2014_preproc(SOFA_obj);
    % Run virtual experiments (demo_reijniers2014)
    num_exp = 100;
    [doa, params] = reijniers2014(template_loc, target, 'num_exp', num_exp);       
    % Calculate performance measures (demo_reijniers2014)
    results.lat_acc = reijniers2014_metrics(doa, 'accL'); % mean lateral error
    results.lat_prec = reijniers2014_metrics(doa, 'precL'); % lateral std
    results.pol_acc = reijniers2014_metrics(doa, 'accP'); % mean polar error
    results.pol_prec = reijniers2014_metrics(doa, 'precP'); % polar std
    results.template_loc = template_loc; % template DTF
    fprintf('\t\tFinished running reijniers2014 (took %0.2f s)...\n',toc);

    %% Run baumgartner2020
    fprintf('\t\tRunning baumgartner2020...\n'); tic
    % Make DTF for median plane directions only
    dtf = getDTF(h(:,indMP,:),fs);
    ndirs = numel(elMP);
    results.ext = nan(ndirs,1);
    template_ext = cell(ndirs,1);
    for j=1:ndirs
        template_ext{j} = makesofa(dtf(:,j,:),fs,azMP(j),elMP(j));
        % Get externalisation values
        results.ext(j) = baumgartner2020(template_ext{j},template_ext{j});
    end
    results.template_ext = template_ext; % template DTFs
    fprintf('\t\tFinished running baumgartner2020 (took %0.2f s)...\n',toc);

    %% Run jelfs2011
    fprintf('\t\tRunning jelfs2011...\n'); tic
    ndirs = numel(azHP);
    srm = nan(ndirs,1);
    target = squeeze(h(:,indFront,:)); % target fixed at front
    for j = 1:ndirs
        interferer = squeeze(h(:,indHP(j),:)); % interferer moves around the HP
        srm(j) = jelfs2011(target,interferer,fs);
    end 
    results.srm = srm;
    fprintf('\t\tFinished running jelfs2011 (took %0.2f s)...\n',toc);

    %% Save results
    fprintf('\t\tSaving results in %s...\n',filename)
    save(filename,'results')
end

mag_ref = results.mag; % to calculate magnitude error
pd_ref = results.pd; % to calculate phase delay error
Lref = results.L_perDirection; % to calculate PSD
Lavg_ref = results.Lavg; % to level-match loudness
calibrationGain = results.calibrationGain; % to calibrate perceptual spectrum
template_loc = results.template_loc; % template localisation model
template_ext = results.template_ext; % template externalisation model
clear results

%% Then, for all test conditions
for N=N_vec % iterate through spatial orders
    fprintf('Processing order %d...\n',N);
    Y = [];
    for i=1:ncond % iterate through conditions
        name = test_conditions{i}.name;
        filename = sprintf('%s/results/ord%0.2d_%s.mat',workdir,N,name);
        if isfile(filename)
            fprintf('\tFound %s. Skipping...\n',filename)
        else
            fprintf('\tGenerating results...\n')

            if isempty(Y)
                Y = AKsh(N, [], az*180/pi, el*180/pi, 'real').';
            end

            %% Load hnm and interpolate to test grid
            hnm = load(sprintf('%s/hnm/ord%0.2d_%s.mat',workdir,N,name)).hnm;                
            isaligned = strcmp(name,'TA');
            hInterp = fromSH(hnm,fs,az,el,isaligned,r);
            HInterp = ffth(hInterp);

            %% Numerical analysis
            % Magnitude and phase delay error 110-point Lebedev grid
            mag = 20*log10(abs(HInterp(:,indLeb,:)));
            pd = (-unwrap(angle(HInterp(:,indLeb,1)))+unwrap(angle(HInterp(:,indLeb,2))))./(2*pi*f)*1e6;
            results.err_mag = mean(abs(mag-mag_ref),2); % avg abs difference across directions
            results.err_pd = mean(abs(pd-pd_ref),2); % avg abs difference across directions
            % Loudness across directions
            [L,wERB] = perceptualSpectrum(HInterp,fs,calibrationGain);
            Lavgd = sum(L.*wERB,1); % ERB-weighted avg loudness per direction
            Lavg = sum(Lavgd(:,indLeb,:).*wLeb.',2); % weighted avg on the Lebedev grid
            L = L - Lavg + Lavg_ref; % level-match with reference
            Lavgd = sum(L.*wERB,1); % recalculate avg loudness per direction
            Ldif = L-Lref; % loudness difference
            PSD = sum(abs(L-Lref).*wERB,1); % avg abs difference per direction (PSD)
            PSDavg = sum(PSD(:,indLeb,:).*wLeb.',2); % avg PSD over the Lebedev grid
            results.L = Lavgd;
            results.PSD = PSD;
            results.PSDavg = PSDavg;
            % Interaural differences for horizontal plane
            results.itd = itdestimator(permute(hInterp(:,indHP,:),[2,3,1]),'fs',fs,...
                'MaxIACCe','lp','upper_cutfreq', 3000,...
                'butterpoly', 10)*1e6;
            results.ild = getILD(hInterp(:,indHP,:),fs);

            %% Reijniers 2014
            fprintf('\t\tRunning reijniers2014...\n'); tic
            % Make DTF and SOFA object for Lebedev grid directions only
            dtf = getDTF(hInterp(:,indLeb,:),fs);
            SOFA_obj = makesofa(dtf,fs,azLeb,elLeb);
            % Preprocessing source information (demo_reijniers2014)
            [~, target] = reijniers2014_preproc(SOFA_obj);
            % Run virtual experiments (demo_reijniers2014)
            num_exp = 100;
            [doa, params] = reijniers2014(template_loc, target, 'num_exp', num_exp);       
            % Calculate performance measures (demo_reijniers2014)
            results.lat_acc = reijniers2014_metrics(doa, 'accL'); % mean lateral error
            results.lat_prec = reijniers2014_metrics(doa, 'precL'); % lateral std
            results.pol_acc = reijniers2014_metrics(doa, 'accP'); % mean polar error
            results.pol_prec = reijniers2014_metrics(doa, 'precP'); % polar std
            fprintf('\t\tFinished running reijniers2014 (took %0.2f s)...\n',toc);

            %% Run baumgartner2020
            fprintf('\t\tRunning baumgartner2020...\n'); tic
            % Make DTF for median plane directions only
            dtf = getDTF(hInterp(:,indMP,:),fs);
            ndirs = numel(elMP);
            results.ext = nan(ndirs,1);
            for j=1:ndirs
                target = makesofa(dtf(:,j,:),fs,azMP(j),elMP(j));
                % Get externalisation values
                results.ext(j) = baumgartner2020(target,template_ext{j});
            end
            fprintf('\t\tFinished running baumgartner2020 (took %0.2f s)...\n',toc);

            %% Run jelfs2011
            fprintf('\t\tRunning jelfs2011...\n'); tic
            ndirs = numel(azHP);
            srm = nan(ndirs,1);
            target = squeeze(hInterp(:,indFront,:)); % target fixed at front
            for j = 1:ndirs
                interferer = squeeze(hInterp(:,indHP(j),:)); % interferer moves around the HP
                srm(j) = jelfs2011(target,interferer,fs);
            end 
            results.srm = srm;
            fprintf('\t\tFinished running jelfs2011 (took %0.2f s)...\n',toc);

            %% Save results
            fprintf('\t\tSaving results in %s...\n',filename)
            save(filename,'results')
        end
    end
end

%% Plot Fig. 1 (SH spectra)
fig1size = [17.78 7];
fig1 = figure('units','centimeters','PaperUnits','centimeters',...
    'PaperSize',fig1size,'Renderer','painters',...
    'pos',[2 2 fig1size(1) fig1size(2)],...
    'paperposition',[0 0 fig1size(1) fig1size(2)]);
% Tight subplot
gap = [.02 .008]; % gap between subplots in norm units (height width)
marg_h = [.18 .15]; % figure height margins in norm units (lower upper)
marg_w = [.05 .08]; % figure width margins in norm units (left right)
[ha, pos] = tight_subplot(1,2,gap,marg_h,marg_w);
% Fig. 1a: original HRTF
axes(ha(1)), plotSHenergy(Hnm_ref(:,:,1),fs); % left only
xlabel('f (Hz)'), ylabel('Order (n)')
axpos = get(ha(1),'pos');
annotation(fig1,'textbox',...
    [axpos(1)+axpos(3)-0.2 axpos(2)+axpos(4)-0.09 0.2 0.09],...
    'String',{'a) Original'},...
    'HorizontalAlignment','right',...
    'FontWeight','bold',...
    'FontSize',7,...
    'FitBoxToText','off',...
    'EdgeColor','none');
set(gca,'fontsize',7)
% Fig. 1b: time-aligned HRTF
axes(ha(2)), plotSHenergy(Hnm_TA_ref(:,:,1),fs);
xlabel('f (Hz)'), ylabel(''), yticklabels({})
axpos = get(ha(2),'pos');
annotation(fig1,'textbox',...
    [axpos(1)+axpos(3)-0.2 axpos(2)+axpos(4)-0.09 0.2 0.09],...
    'String',{'b) Time-aligned'},...
    'HorizontalAlignment','right',...
    'FontWeight','bold',...
    'FontSize',7,...
    'FitBoxToText','off',...
    'EdgeColor','none');
set(gca,'fontsize',7)
sgtitle('HRTF energy in the SH domain, with and without alignment')

%% Plot Fig. 2 (mag/phase errors per spatial order)
% Load data from file
names = {
    'ord01_Trunc'
    'ord05_Trunc'
    'ord10_Trunc'
    'ord20_Trunc'
    'ord30_Trunc'
    'ord40_Trunc'
};
labels = {
    'N=1'
    'N=5'
    'N=10'
    'N=20'
    'N=30'
    'N=40'  
};
n = numel(names);
err_mag = zeros(nfreqs,n);
err_pd = zeros(nfreqs,n);
for i=1:n
    name = names{i};
    filename = sprintf('%s/results/%s.mat',workdir,name);
    load(filename,'results');
    err_mag(:,i) = results.err_mag(:,1,1); % left ear only
    err_pd(:,i) = results.err_pd;
end

% Tight subplot
fig2 = figure('pos',[56.6000 98.6000 323.2000 319.6000]);
gap = [.02 .008]; % gap between subplots in norm units (height width)
marg_h = [.09 .1]; % [.09 .1] figure height margins in norm units (lower upper)
marg_w = [.1 .03]; % figure width margins in norm units (left right)
[ha, ~] = tight_subplot(2,1,gap,marg_h,marg_w);

% Top plot: magnitude error
colors = parula(n+1);
lsvec = {'-'};%,'-.'};
lwvec = [0.5,0.5];
mvec = {'^','v','x','s','o','d','p','h'};
ms = 3; % marker size
mi = int32(logspace(log10(1),log10(1025),10));
    
axes(ha(1))
for i=1:n
    ls = lsvec{mod(i-1,numel(lsvec))+1};
    lw = lwvec(mod(i-1,numel(lwvec))+1);
    m = mvec{i};
    semilogx(f,err_mag(:,i),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi); hold on
end

semilogx([f(2) 20000], [1 1],'k:')
legend(labels,'location','west')

ylim([0 40])
xlabel(''), xticklabels({}), xlim([f(2) 20000])
ylabel('Error (dB)'), grid on
axpos = get(ha(1),'pos');
annotation(fig2,'textbox',...
    [axpos(1) axpos(2)+axpos(4)-0.09 0.5 0.09],...
    'String',{'a. Magnitude error'},...
    'HorizontalAlignment','left',...
    'FontWeight','bold',...
    'FontSize',7,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Bottom plot: phase error
axes(ha(2))
for i=1:n
    ls = lsvec{mod(i-1,numel(lsvec))+1};
    lw = lwvec(mod(i-1,numel(lwvec))+1);
    m = mvec{i};
    semilogx(f,err_pd(:,i),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi); hold on
end

semilogx([f(2) 20000], [20 20],'k:')
ylim([0 300])
xlim([f(2) 20000]), grid on, ylabel('Error (\mus)')
xticks([100,1000,10000,20000])
xticklabels(["100","1k","10k","20k"])    
xlabel('f (Hz)')
axpos = get(ha(2),'pos');
annotation(fig2,'textbox',...
    [axpos(1) axpos(2)+axpos(4)-0.09 0.5 0.09],...
    'String',{'b. Phase delay error'},...
    'HorizontalAlignment','left',...
    'FontWeight','bold',...
    'FontSize',7,...
    'FitBoxToText','off',...
    'EdgeColor','none');

set(ha,'fontsize',7)
set(gcf,'units','centimeters','Renderer','painters')
figpos = get(gcf,'position');
set(gcf,'PaperSize',[figpos(3) figpos(4)],'paperposition',[0 0 figpos(3) figpos(4)])
sgtitle('Interpolation errors (non-aligned HRTF)')

%% Plot Fig. 3 (mag/phase errors per spatial order, time-aligned)
% Load data from file
names = {
    'ord01_TA'
    'ord05_TA'
    'ord10_TA'
    'ord20_TA'
    'ord30_TA'
    'ord40_TA'
};
labels = {
    'N=1'
    'N=5'
    'N=10'
    'N=20'
    'N=30'
    'N=40'  
};
n = numel(names);
err_mag = zeros(nfreqs,n);
err_pd = zeros(nfreqs,n);
for i=1:n
    name = names{i};
    filename = sprintf('%s/results/%s.mat',workdir,name);
    load(filename,'results');
    err_mag(:,i) = results.err_mag(:,1,1); % left ear only
    err_pd(:,i) = results.err_pd;
end

% Tight subplot
fig3 = figure('pos',[56.6000 98.6000 323.2000 319.6000]);
gap = [.02 .008]; % gap between subplots in norm units (height width)
marg_h = [.09 .1]; % [.09 .1] figure height margins in norm units (lower upper)
marg_w = [.1 .03]; % figure width margins in norm units (left right)
[ha, ~] = tight_subplot(2,1,gap,marg_h,marg_w);

% Top plot: magnitude error
colors = parula(n+1);
lsvec = {'-'};%,'-.'};
lwvec = [0.5,0.5];
mvec = {'^','v','x','s','o','d','p','h'};
ms = 3; % marker size
mi = int32(logspace(log10(1),log10(1025),10));
    
axes(ha(1))
for i=1:n
    ls = lsvec{mod(i-1,numel(lsvec))+1};
    lw = lwvec(mod(i-1,numel(lwvec))+1);
    m = mvec{i};
    semilogx(f,err_mag(:,i),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi); hold on
end

semilogx([f(2) 20000], [1 1],'k:')
legend(labels,'location','west')

ylim([0 40])
xlabel(''), xticklabels({}), xlim([f(2) 20000])
ylabel('Error (dB)'), grid on
axpos = get(ha(1),'pos');
annotation(fig2,'textbox',...
    [axpos(1) axpos(2)+axpos(4)-0.09 0.5 0.09],...
    'String',{'a. Magnitude error'},...
    'HorizontalAlignment','left',...
    'FontWeight','bold',...
    'FontSize',7,...
    'FitBoxToText','off',...
    'EdgeColor','none');

% Bottom plot: phase error
axes(ha(2))
for i=1:n
    ls = lsvec{mod(i-1,numel(lsvec))+1};
    lw = lwvec(mod(i-1,numel(lwvec))+1);
    m = mvec{i};
    semilogx(f,err_pd(:,i),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi); hold on
end

semilogx([f(2) 20000], [20 20],'k:')
ylim([0 300])
xlim([f(2) 20000]), grid on, ylabel('Error (\mus)')
xticks([100,1000,10000,20000])
xticklabels(["100","1k","10k","20k"])    
xlabel('f (Hz)')
axpos = get(ha(2),'pos');
annotation(fig2,'textbox',...
    [axpos(1) axpos(2)+axpos(4)-0.09 0.5 0.09],...
    'String',{'b. Phase delay error'},...
    'HorizontalAlignment','left',...
    'FontWeight','bold',...
    'FontSize',7,...
    'FitBoxToText','off',...
    'EdgeColor','none');

set(ha,'fontsize',7)
set(gcf,'units','centimeters','Renderer','painters')
figpos = get(gcf,'position');
set(gcf,'PaperSize',[figpos(3) figpos(4)],'paperposition',[0 0 figpos(3) figpos(4)])
sgtitle('Interpolation errors (time-aligned HRTF)')

%% Plot Fig. 3 (ITD/ILD per spatial order)
% Load data from file
names = {
    'ord01_Trunc'
    'ord05_Trunc'
    'ord10_Trunc'
    'ord20_Trunc'
    'ord30_Trunc'
    'ord40_Trunc'
    'ref'
};
labels = {
    'N=1'
    'N=5'
    'N=10'
    'N=20'
    'N=30'
    'N=40'
    'Reference'
};
n = numel(names);
ndirs = numel(azHP);
itd = zeros(ndirs,n);
ild = zeros(ndirs,n);
for i=1:n
    name = names{i};
    filename = sprintf('%s/results/%s.mat',workdir,name);
    load(filename,'results');
    itd(:,i) = results.itd;
    ild(:,i) = results.ild;
end

% Tight subplot  
fig4size = [21 13];
fig4 = figure('units','centimeters','pos',[2 2 fig4size(1) fig4size(2)],...
    'Renderer','painters','PaperSize',[fig4size(1) fig4size(2)],...
    'paperposition',[0 0 fig4size(1) fig4size(2)]);
gap = [.04 .06]; % gap between subplots in norm units (height width)
marg_h = [0.05 .15]; % [.09 .1] figure height margins in norm units (lower upper)
marg_w = [.07 .03]; % figure width margins in norm units (left right)
[ha, ~] = tight_subplot(2,2,gap,marg_h,marg_w); 
 
colors = parula(n);
lsvec = {'-'};%,'-.'};
lwvec = [0.5,0.5];
mvec = {'^','v','x','s','o','d','p','h'};
ms = 3; % marker size
step_big = round(numel(azHP)/4);
step_small = round(step_big/n);
    
% ITD
axes(ha(1))
for i=1:n % multiply by 0.001 to have it in ms (cleaner plot)
    ls = lsvec{mod(i-1,numel(lsvec))+1};
    lw = lwvec(mod(i-1,numel(lwvec))+1);
    m = mvec{i};
    mi = [((i-1)*step_small+1):step_big:numel(azHP)-1]; % marker indices
    color = colors(i,:);
    if i==n % reference
        ls = ':';
        lw = 0.5;
        m = 'none';
        color = [0 0 0];
    end
    polarplot(azHP,abs(itd(:,i))*0.001,'Color',color,'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi), hold on

end
set(gca,'ThetaZeroLocation','top','ThetaLim',[0 360],'fontsize',7)
set(gca,'ThetaTick',[0:45:360],'ThetaTickLabel',["0","45","90","135","","-135","-90","-45"],'RAxisLocation',-90)
title('ITD (ms)')

% ILD
axes(ha(2))
for i=1:n
    ls = lsvec{mod(i-1,numel(lsvec))+1};
    lw = lwvec(mod(i-1,numel(lwvec))+1);
    m = mvec{i};
    mi = [((i-1)*step_small+1):step_big:numel(azHP)-1]; % marker indices
    color = colors(i,:);
    if i==n % reference
        ls = ':';
        lw = 0.5;
        m = 'none';
        color = [0 0 0];
    end
    polarplot(azHP,abs(ild(:,i)),'Color',color,'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi), hold on
end
set(gca,'ThetaZeroLocation','top','ThetaLim',[0 360],'fontsize',7)
set(gca,'ThetaTick',[0:45:360],'ThetaTickLabel',["0","45","90","135","","-135","-90","-45"],'RAxisLocation',-90)
title('ILD (dB)')

legend(labels,'position',[0.4646    0.6219    0.1135    0.1883]);

itderr = abs(itd-itd(:,end));
ilderr = abs(ild-ild(:,end));

axes(ha(3)), violinplot(itderr(:,1:end-1),[],'ShowData',false,'BoxWidth',0.03); grid on
set(gca,'YTickLabelMode','auto')
hold on, plot([0 n],[itdJND itdJND],'k:')
xticklabels(labels(1:end-1))
ylabel('Abs. ITD error (ms)')
set(gca,'fontsize',6)
ylim([0 700])

axes(ha(4)), violinplot(ilderr(:,1:end-1),[],'ShowData',false,'BoxWidth',0.03); grid on
set(gca,'YTickLabelMode','auto')
hold on, plot([0 n],[ildJND ildJND],'k:')
xticklabels(labels(1:end-1))
ylabel('Abs. ILD error (dB)')
set(gca,'fontsize',6)
ylim([0 10])

sgtitle('Interaural differences (non-aligned HRTF)')

%% Plot Fig. 5 (ITD/ILD per spatial order)
% Load data from file
names = {
    'ord01_TA'
    'ord05_TA'
    'ord10_TA'
    'ord20_TA'
    'ord30_TA'
    'ord40_TA'
    'ref'
};
labels = {
    'N=1'
    'N=5'
    'N=10'
    'N=20'
    'N=30'
    'N=40'
    'Reference'
};
n = numel(names);
ndirs = numel(azHP);
itd = zeros(ndirs,n);
ild = zeros(ndirs,n);
for i=1:n
    name = names{i};
    filename = sprintf('%s/results/%s.mat',workdir,name);
    load(filename,'results');
    itd(:,i) = results.itd;
    ild(:,i) = results.ild;
end

% Tight subplot  
fig5 = figure('units','centimeters','pos',[2 2 fig4size(1) fig4size(2)],...
    'Renderer','painters','PaperSize',[fig4size(1) fig4size(2)],...
    'paperposition',[0 0 fig4size(1) fig4size(2)]);
gap = [.04 .06]; % gap between subplots in norm units (height width)
marg_h = [0.05 .15]; % [.09 .1] figure height margins in norm units (lower upper)
marg_w = [.07 .03]; % figure width margins in norm units (left right)
[ha, ~] = tight_subplot(2,2,gap,marg_h,marg_w); 
 
colors = parula(n);
lsvec = {'-'};%,'-.'};
lwvec = [0.5,0.5];
mvec = {'^','v','x','s','o','d','p','h'};
ms = 3; % marker size
step_big = round(numel(azHP)/4);
step_small = round(step_big/n);
    
% ITD
axes(ha(1))
for i=1:n % multiply by 0.001 to have it in ms (cleaner plot)
    ls = lsvec{mod(i-1,numel(lsvec))+1};
    lw = lwvec(mod(i-1,numel(lwvec))+1);
    m = mvec{i};
    mi = [((i-1)*step_small+1):step_big:numel(azHP)-1]; % marker indices
    color = colors(i,:);
    if i==n % reference
        ls = ':';
        lw = 0.5;
        m = 'none';
        color = [0 0 0];
    end
    polarplot(azHP,abs(itd(:,i))*0.001,'Color',color,'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi), hold on

end
set(gca,'ThetaZeroLocation','top','ThetaLim',[0 360],'fontsize',7)
set(gca,'ThetaTick',[0:45:360],'ThetaTickLabel',["0","45","90","135","","-135","-90","-45"],'RAxisLocation',-90)
title('ITD (ms)')

% ILD
axes(ha(2))
for i=1:n
    ls = lsvec{mod(i-1,numel(lsvec))+1};
    lw = lwvec(mod(i-1,numel(lwvec))+1);
    m = mvec{i};
    mi = [((i-1)*step_small+1):step_big:numel(azHP)-1]; % marker indices
    color = colors(i,:);
    if i==n % reference
        ls = ':';
        lw = 0.5;
        m = 'none';
        color = [0 0 0];
    end
    polarplot(azHP,abs(ild(:,i)),'Color',color,'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi), hold on
end
set(gca,'ThetaZeroLocation','top','ThetaLim',[0 360],'fontsize',7)
set(gca,'ThetaTick',[0:45:360],'ThetaTickLabel',["0","45","90","135","","-135","-90","-45"],'RAxisLocation',-90)
title('ILD (dB)')

legend(labels,'position',[0.4646    0.6219    0.1135    0.1883]);

itderr = abs(itd-itd(:,end));
ilderr = abs(ild-ild(:,end));

axes(ha(3)), violinplot(itderr(:,1:end-1),[],'ShowData',false,'BoxWidth',0.03); grid on
set(gca,'YTickLabelMode','auto')
hold on, plot([0 n],[itdJND itdJND],'k:')
xticklabels(labels(1:end-1))
ylabel('Abs. ITD error (ms)')
set(gca,'fontsize',6)
ylim([0 700])

axes(ha(4)), violinplot(ilderr(:,1:end-1),[],'ShowData',false,'BoxWidth',0.03); grid on
set(gca,'YTickLabelMode','auto')
hold on, plot([0 n],[ildJND ildJND],'k:')
xticklabels(labels(1:end-1))
ylabel('Abs. ILD error (dB)')
set(gca,'fontsize',6)
ylim([0 10])

sgtitle('Interaural differences (time-aligned HRTF)')

%% Plot Fig. 6 (loudness map)
% Load data from file
names = {
    'ref'
    'ord05_Trunc'
    'ord05_TA'
};
labels = {
    'Reference'
    'Trunc'
    'TA'
};
n = numel(names);
ndirs = numel(az);
L = cell(n,1);
PSD = zeros(ndirs,n);
PSDavg = zeros(n,1);
for i=1:n
    name = names{i};
    filename = sprintf('%s/results/%s.mat',workdir,name);
    load(filename,'results');
    L{i} = results.L;
    if isfield(results,'PSD')
        PSD(:,i) = results.PSD(:,:,1).'; % left ear only
        PSDavg(i) = results.PSDavg(:,:,1); % left ear only
    end
end

% Tight subplot
fig6size = [16.4253 8];
fig6 = figure('units','centimeters','pos',[2 2 fig6size(1) fig6size(2)],...
    'Renderer','painters','PaperSize',[fig6size(1) fig6size(2)],...
    'paperposition',[0 0 fig6size(1) fig6size(2)]);
gap = [.07 .008]; % gap between subplots in norm units (height width)
marg_h = [.5 .15]; % [.09 .1] figure height margins in norm units (lower upper)
marg_w = [.06 .1]; % figure width margins in norm units (left right)
[ha, ~] = tight_subplot(1,3,gap,marg_h,marg_w);
clims = [1 5]; % change as needed

colormap parula
for i=1:n
    axes(ha(i)) 
    plotSph(az,el,L{i}(:,:,1)) % left ear only
    if i==1 || i==6
        yticks(-60:30:60)
        ylabel('Elevation (deg)')
    else
        yticklabels({})
    end
    xlabel('Azimuth (deg)')
    if PSDavg(i)>0
        title(sprintf('%s [%0.2f]',labels{i},PSDavg(i)));
    else
        title(labels{i}) % reference
    end

    set(gca,'fontsize',7,'XDir','reverse')
    grid(gca,'off')
    caxis(clims)
    if i==n
        c = colorbar; c.Label.String = 'Loudness (sones)';
        c.Position = [0.9174    0.4904    0.0111    0.3804];
    end
end 

ax = axes('innerposition',[0.4 0.05 0.2 0.3]);
violinplot(PSD(:,2:end),[],'ShowData',false,'BoxWidth',0.03); grid on
% boxplot(PSD(:,2:end)), grid on
% hold on, plot(mean(PSD(:,2:end)),'dg')
xticklabels(labels(2:end))
ylabel('PSD (sones)')
set(gca,'fontsize',7)

sgtitle('Loudness maps and perceptual spectral difference (PSD)')

%% Plot Fig. 7 (models' outputs per spatial order)
% Load data from file
names = {
    'Trunc'
    'TA'
    'ref'
};
labels = {
    'Trunc'
    'TA'
    'Reference'
};
n = numel(names);
m = numel(N_vec);
PSD = nan(n,m);
lat_prec = nan(n,m);
pol_prec = nan(n,m);
ext = nan(n,m);
srm = nan(n,m);
for i=1:n
    name = names{i};
    if strcmp(name,'ref')
        filename = sprintf('%s/results/ref.mat',workdir);
        load(filename,'results');
        PSD(i,:) = 0;
        lat_prec(i,:) = results.lat_prec;
        pol_prec(i,:) = results.pol_prec;
        ext(i,:) = mean(results.ext); 
        srm(i,:) = mean(results.srm);
    else
        for j=1:m
            N=N_vec(j);
            filename = sprintf('%s/results/ord%0.2d_%s.mat',workdir,N,name);
            load(filename,'results');
            PSD(i,j) = results.PSDavg(:,:,1); % left ear only
            lat_prec(i,j) = results.lat_prec;
            pol_prec(i,j) = results.pol_prec;
            ext(i,j) = mean(results.ext); 
            srm(i,j) = mean(results.srm);
        end
    end
end
% Tight subplot
fig7size = [18.3938 9.7790];
fig7 = figure('units','centimeters','pos',[2 2 fig7size(1) fig7size(2)],...
    'Renderer','painters','PaperSize',[fig7size(1) fig7size(2)],...
    'paperposition',[0 0 fig7size(1) fig7size(2)]);
gap = [.06 .05]; % gap between subplots in norm units (height width)
marg_h = [.1 .1]; % [.09 .1] figure height margins in norm units (lower upper)
marg_w = [.07 .02]; % figure width margins in norm units (left right)
[ha, ~] = tight_subplot(2,3,gap,marg_h,marg_w); 
    
colors = parula(n); 
lsvec = {'-'};%,'-.'};
mvec = {'^','v','x','s','o','d','p','h'};
lwvec = [0.5,0.5];
ms = 3; % marker size
mi = 1:n; % [1:5:44]; % marker indices
    
colors(end,:) = [0 0 0]; % last one is reference
  
for i=1:n
    if i<n
        ls = lsvec{mod(i-1,numel(lsvec))+1};
        lw = lwvec(mod(i-1,numel(lwvec))+1);
        m = mvec{i};
    else
        ls = ':'; % reference
        lw = 0.5;
        m = 'None';
    end
    plot(ha(1),N_vec,PSD(i,:),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi),hold(ha(1),'on')
    plot(ha(2),N_vec,lat_prec(i,:),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi),hold(ha(2),'on')
    plot(ha(3),N_vec,pol_prec(i,:),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi),hold(ha(3),'on')
    plot(ha(4),N_vec,ext(i,:),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi),hold(ha(4),'on')
    plot(ha(5),N_vec,srm(i,:),'Color',colors(i,:),'LineWidth',lw,'LineStyle',ls,'Marker',m,'MarkerSize',ms,'MarkerIndices',mi),hold(ha(5),'on')
end
legh = legend(ha(1),labels,'location','best');

ylabel(ha(1),'PSD (sones)')%, ylim(ha(1),[min(PSD(:),max(PSD(:)])
ylabel(ha(2),'Lateral precision (deg)')
ylabel(ha(3),'Polar precision (deg)')
ylabel(ha(4),'Externalisation')%, ylim(ha(4),[0.3 1.01])
ylabel(ha(5),'SRM (dB)') % , ylim(ha(4),[0.3 1.01])
for i=1:5
    grid(ha(i),'on')
    xlim(ha(i),[1 44])
    xticks(ha(i),[1,5:5:44])
    if i<=3
%             xticklabels(ha(i),{})
    else
        xlabel(ha(i),'Spatial order (N)')
    end
    set(ha(i),'fontsize',7)
end
set(ha(6),'visible','off')

sgtitle('Perceptual models'' output (per order)')

