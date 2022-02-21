function SOFA_obj = hrtf2sofa(h,fs,az,el,r)
% Take HRIR data matrix and other paramters and make a SOFA object. Assume
% "normal" conditions (free field, listener at the origin looking at the
% front, etc).
%
% INPUT:
%   h = HRIRs (time x dirs x ears)
%   fs = sampling frequency in Hz
%   az = azimuth in rad (dirs x 1)
%   el = colatitude in rad (dirs x 1)
%   r = head radius in m (def=0.085)
%
% OUTPUT:
%   SOFA_obj = SOFA object, e.g. obtained from SOFAload()
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% February 2021

if ~exist('r','var')
    r = 0.085;
end

SOFA_obj = SOFAgetConventions('SimpleFreeFieldHRIR');
SOFA_obj.Data.IR = permute(h,[2,3,1]); % dirs x ears x time
SOFA_obj.Data.SamplingRate = fs;
SOFA_obj.Data.SamplingRate_Units = 'hertz';
SOFA_obj.SourcePosition = [az(:)*180/pi,90-el(:)*180/pi,ones(size(az(:)))]; % HRTF directions
SOFA_obj.SourcePosition_Type = 'spherical';
SOFA_obj.SourcePosition_Units = 'degree,degree,metre';
SOFA_obj.ListenerPosition = [0 0 0]; % listener at the origin
SOFA_obj.ListenerPosition_Type = 'cartesian';
SOFA_obj.ListenerPosition_Units = 'metre';
SOFA_obj.ListenerView = [1 0 0]; % listener looking front
SOFA_obj.ListenerView_Type = 'cartesian';
SOFA_obj.ListenerView_Units = 'metre';
SOFA_obj.ReceiverPosition = [0 r 0; 0 -r 0];
SOFA_obj.ReceiverPosition_Type = 'cartesian';
SOFA_obj.ReceiverPosition_Units = 'metre';
% SOFA_obj.EmitterPosition = [0,0,0]; % not used; this is the default value
% SOFA_obj.EmitterPosition_Type = 'cartesian';
% SOFA_obj.EmitterPosition_Units = 'metre';
SOFA_obj.API.M = size(SOFA_obj.Data.IR,1);
SOFA_obj.API.R = size(SOFA_obj.Data.IR,2);
SOFA_obj.API.N = size(SOFA_obj.Data.IR,3);
SOFA_obj.API.Dimensions.Data.IR = 'MRN'; % time-measurements-receivers
SOFA_obj.GLOBAL_Comment = 'Generated with ''hrtf2sofa()'' from an HRIR';
SOFA_obj.GLOBAL_Title = 'HRTF';