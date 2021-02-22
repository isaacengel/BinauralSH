function [Hnm,G1,G2] = applyCovConst(Hnm,Hnm_ref,debug_plots,fs)
% Apply diffuse field and covariance constraint as described in [1] sec.
% 4.11.3, which refers to [2] and [3].
%
% Important: the method works better if the HRIRs are aligned at t=0. This
% is taken care of in the function 'toSH'.
%
% INPUT:
%   Hnm = HRTF's SH coefficients (nfreqs x (N+1)^2 x 2 ears)
%   Hnm_ref = reference (high-order) HRTF's SH coefficients
%   debug_plots = whether to show some debug plots (def=false)
%   fs = sampling rate in Hz
%
% OUTPUT:
%   Hnm = equalised HRTF's SH coefficients
%   G1 = global EQ filter (nfreqs x 2 ears)
%   G2 = left/right mix (nfreqs x 2 ears)
%
% REFERENCES:
%   [1] Zotter, Franz, and Matthias Frank. Ambisonics: A practical 3D audio
%       theory for recording, studio production, sound reinforcement, and
%       virtual reality. Springer Nature, 2019.
%   [2] Zaunschirm, Markus, Christian Schörkhuber, and Robert Höldrich.
%       "Binaural rendering of Ambisonic signals by head-related impulse
%       response time alignment and a diffuseness constraint." The Journal
%       of the Acoustical Society of America 143.6 (2018): 3616-3627.
%   [3] Vilkamo, Juha, Tom Bäckström, and Achim Kuntz. "Optimized
%       covariance domain framework for time–frequency processing of
%       spatial audio." Journal of the Audio Engineering Society 61.6
%       (2013): 403-411.
%
% AUTHOR: Isaac Engel (isaac.engel@imperial.ac.uk)
% February 2021

%% Check inputs
if ~exist('debug_plots','var')
    debug_plots = 0;
end

%% Some parameters
nfreqs = size(Hnm,1); 
nsh = size(Hnm,2);
nears = size(Hnm,3);
N = sqrt(nsh)-1;
G1 = zeros(ones,nears); % diffuse EQ (nfreqs x nears)
G2 = zeros(nfreqs,nears); % mix-channel EQ (nfreqs x nears)
if N<1
    warning('Hnm order was 0. Covariance constraint not applied...')
    return
end
if nears < 2
    warning('Hnm was one ear only. Covariance constraint not applied...')
end

%% Covariance constraint ([2] sec.4.11.3)
fprintf('Applying covariance constraint')
Hnm_match = zeros(size(Hnm));
count = 0; % for text output
if debug_plots % variables for debugging
    df_ref = zeros(nfreqs,2);
    df_1 = zeros(nfreqs,2);
    df_2 = zeros(nfreqs,2);
    ic_ref = zeros(nfreqs,1);
    ic_1 = zeros(nfreqs,1);
    ic_2 = zeros(nfreqs,1);
    EQprod = zeros(nfreqs,nears);
    EQsum = zeros(nfreqs,nsh,nears);
end
for fi=1:nfreqs
    Hnmref_k = squeeze(Hnm_ref(fi,:,:)); % reference HRTF at fi [ndirs x 2]
    Hnm_k = squeeze(Hnm(fi,:,:)); % interpolated HRTF at fi [ndirs x 2]
    Rref = Hnmref_k'*Hnmref_k/(4*pi);%size(Hnmref_k,1); % reference HRTF covariance matrix
    R = Hnm_k'*Hnm_k/(4*pi);%size(Hnm_k,1); % interp. HRTF covariance matrix
    Rref(1,1) = real(Rref(1,1)); % force diagonals to be real
    Rref(2,2) = real(Rref(2,2)); % Imag part is due to numerical error
    R(1,1) = real(R(1,1)); 
    R(2,2) = real(R(2,2)); 
    Xref = chol(Rref); % Xref'*Xref == Href_k'*Href_k
    X = chol(R); % X'*X == H_k'*H_k
    [U,S,V] = svd(X'*Xref); % [1] eq.A.72 and above
    if ~isreal(S) || any(S(:)<0)
        warning('The SVD decomposition yielded negative or complex values')
    end
    M = X\(V*U'*Xref);
    Hnm_match(fi,:,:) = squeeze(Hnm(fi,:,:))*M; % [1] eq.4.63 and eq.A.72
    G1(fi,:) = [M(1,1),M(2,2)]; % diffuse EQ at fi (L,R)
    G2(fi,:) = [M(2,1),M(1,2)]; % mix-channel EQ at fi (L,R)
    
    % Update text output
    count = count+1;
    if count >= nfreqs/10
        count = 0;
        fprintf('.');
    end
    if debug_plots
        % Calculate covariance matrix for the corrected low-order HRTF
        Hmatch_k = squeeze(Hnm_match(fi,:,:));
        Rm = Hmatch_k'*Hmatch_k/(4*pi);
        Rm(1,1) = real(Rm(1,1)); 
        Rm(2,2) = real(Rm(2,2)); 
        % Save data to plot later
        df_ref(fi,:) = [Rref(1,1),Rref(2,2)];
        df_1(fi,:) = [Rm(1,1),Rm(2,2)];
        df_2(fi,:) = [R(1,1),R(2,2)];
        ic_ref(fi) = abs(Rref(1,2))/sqrt(Rref(1,1)*Rref(2,2));
        ic_1(fi) = abs(Rm(1,2))/sqrt(Rm(1,1)*Rm(2,2));
        ic_2(fi) = abs(R(1,2))/sqrt(R(1,1)*R(2,2));
        EQprod(fi,:) = [M(1,1),M(2,2)];
        EQsum(fi,:,:) = cat(3,M(2,1)*squeeze(Hnm(fi,:,2)),M(1,2)*squeeze(Hnm(fi,:,1)));
    end
end
fprintf('\n')

% Hnm_match2(:,:,1) = Hnm(:,:,1).*G1(:,1) + Hnm(:,:,2).*G2(:,1);
% Hnm_match2(:,:,2) = Hnm(:,:,2).*G1(:,2) + Hnm(:,:,1).*G2(:,2);
% 
% Hnm_match3 = Hnm.*permute(G1,[1,3,2]);
% Hnm_match3(:,:,1) = Hnm_match3(:,:,1) + Hnm(:,:,2).*G2(:,1);
% Hnm_match3(:,:,2) = Hnm_match3(:,:,2) + Hnm(:,:,1).*G2(:,2);


Hnm = Hnm_match;

%% Plots
% Test diffuse field and IC compared to ref
% This is trying to reproduce Fig. 7 from [2]
if debug_plots
    f = linspace(0,fs/2,nfreqs);
    figure
    subplot(2,2,1)
    semilogx(f,db(abs(df_ref(:,1))),'k--','Linewidth',1.5), hold on
    semilogx(f,db(abs(df_1(:,1))))
    semilogx(f,db(abs(df_2(:,1))))
    xlim([20 20000]), grid on, xlabel('f (Hz)'), ylabel('dB')
    title('diffuse field energy left')
    subplot(2,2,2)
    semilogx(f,db(abs(df_ref(:,2))),'k--','Linewidth',1.5), hold on
    semilogx(f,db(abs(df_1(:,2))))
    semilogx(f,db(abs(df_2(:,2))))
    xlim([20 20000]), grid on, xlabel('f (Hz)'), ylabel('dB')
    title('diffuse field energy right')
    subplot(2,2,3)
    semilogx(f,ic_ref,'k--','Linewidth',1.5), hold on
    semilogx(f,ic_1)
    semilogx(f,ic_2)
    xlim([20 20000]), ylim([0 1]), grid on, xlabel('f (Hz)'), ylabel('IACC')
    legend({'Ref','matched conv const','not matched'},'position',[0.5972,0.2236,0.2632,0.1090])
    title('Interaural coherence')
    sgtitle(sprintf('N=%d',N))
end
