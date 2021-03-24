function [h,fs,az,el,r] = sofa2hrtf(SOFA_obj)
% Return HRIR data matrix and other paramters from a SOFA object. Assume
% "normal" conditions (free field, listener at the origin looking at the
% front, etc).
%
% INPUT:
%   SOFA_obj = SOFA object, e.g. obtained from SOFAload()
%
% OUTPUT:
%   h = HRIRs (time x dirs x ears)
%   fs = sampling frequency in Hz
%   az = azimuth in rad (dirs x 1)
%   el = colatitude in rad (dirs x 1)
%   r = head radius in m (def=0.085)
%
% EXTERNAL DEPENDENCIES:
%   SOFA API for Matlab (github.com/sofacoustics/API_MO)
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021

%% Check struct fields
if ~strcmp(SOFA_obj.GLOBAL_SOFAConventions,'SimpleFreeFieldHRIR')
    error('SOFA convention must be SimpleFreeFieldHRIR')
end
if ~strcmp(SOFA_obj.SourcePosition_Type,'spherical')
    error('Source position must be in spherical coordinates')
end
if length(SOFA_obj.API.Dimensions.Data.IR)~=3
    error('IR must have 3 dimensions (direction, ear, time)')
end

%% Assign data to output
h = SOFA_obj.Data.IR;
fs = SOFA_obj.Data.SamplingRate;
az = SOFA_obj.SourcePosition(:,1)*pi/180; % azimuth in rad
el = pi/2-SOFA_obj.SourcePosition(:,2)*pi/180; % colatitude in rad
r = 0.0875; % using this as default for now

%% Rearrange IR dimensions
dimIn = SOFA_obj.API.Dimensions.Data.IR; % dimensions, e.g. 'MRN'
dimOut = 'NMR'; % this is the order we want for the output
dimIndx = nan(1,3);
for i=1:3 % this is the order used in the output
    dimIndx(i) = find(dimIn==dimOut(i));
end
h = permute(h,dimIndx);
