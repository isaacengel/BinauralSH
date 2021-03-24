function h = fromSH(hnm,fs,az,el,isaligned,r)
% Peform SH interpolation to obtain HRIRs at discrete directions from their
% SH representation.
%
% SIMPLE USAGE EXAMPLE:
%   h = fromSH(hnm,az,el);
%
% INPUT:
%   hnm = SH coefficients of HRIR set (irlen x (N+1)^2 x 2 ears)
%   fs = sampling frequency in Hz
%   az = HRIR azimuth (ndirs x 1) in rad (0=front, pi/2=left)
%   el = HRIR elevation (ndirs x 1) in rad (0=top, pi/2=front)
%   isaligned = whether the HRIRs have undergone time-alignment (see
%       toSH_TA; default=0)
%   r = head radius in m used for the alignment (default=0.0875)
%
% OUTPUT:
%   h = HRIRs in matrix format (irlen x ndirs x 2 ears)
%
% EXTERNAL DEPENDENCIES:
%   AKtools (www.ak.tu-berlin.de/aktools)
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% March 2021

%% Check inputs
if ~exist('isaligned','var')
    isaligned = 0;
end
if ~exist('r','var')
    r = 0.0875;
end

%% Get relevant parameters
[irlen,nsh,nears] = size(hnm);
nfft = 2^nextpow2(irlen);
nfreqs = nfft/2+1;
f = linspace(0,fs/2,nfreqs).'; % frequency vector
N = sqrt(nsh)-1;
c = 343; % c=343m/s
k = 2*pi*f.'/c; % wave number
kr = k*r;

%% Process
Y = AKsh(N, [], az*180/pi, el*180/pi, 'real').'; % SH coefficients
Hnm = ffth(hnm,nfft,1); % to frequency domain
H = pagemtimes(Hnm,Y); % interpolate
if isaligned % re-align if specified
    p = earAlign(kr,az,el);
    H = H.*exp(1i*p);
end
h = iffth(H); % to time domain

