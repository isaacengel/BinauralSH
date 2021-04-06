function p = earAlign(kr,az,el,earAz,earEl)
% Return phase corrections for ear alignment [1], such that:
%   H_aligned = H.*exp(-1i*p);
%   H = H_aligned.*exp(1i*p);
%
% INPUT:
%   kr = kr vector = 2*pi*f*r/c (nfreqs x 1)
%   az = HRIR azimuth (ndirs x 1) in rad
%   el = HRIR elevation (ndirs x 1) in rad (0=top, pi/2=front)
%   earAz = left/right ear azimuth (1 x 2) in rad (def=[pi/2, 3*pi/2])
%   earEl = left/right ear elevation (1 x 2) in rad (def = [pi/2, pi/2])
%   
% OUTPUT:
%   p = phase correction (nfreqs x ndirs x 2 ears)
%
% REFERENCES:
%   [1] Ben-Hur, Zamir, et al. "Efficient Representation and Sparse
%       Sampling of Head-Related Transfer Functions Using Phase-Correction
%       Based on Ear Alignment." IEEE/ACM Transactions on Audio, Speech,
%       and Language Processing 27.12 (2019): 2249-2262.
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% March 2021

%% Prepare parameters
if ~exist('earAz','var') || isempty(earAz)
    earAz = [pi/2, 3*pi/2];
end
if ~exist('earEl','var') || isempty(earEl)
    earEl = [pi/2, pi/2];
end
az = az(:).'; % force row vector
el = el(:).';
kr = kr(:); % force column vector

%% Process
costhetal = cos(el)*cos(earEl(1))+sin(el)*sin(earEl(1)).*cos(az-earAz(1));
costhetar = cos(el)*cos(earEl(2))+sin(el)*sin(earEl(2)).*cos(az-earAz(2));
p = cat(3,kr*costhetal,kr*costhetar); % phase correction
