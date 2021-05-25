function binauralSH_stop()
%binauralSH_start Stop the Binaural SH Toolbox
%
%   Usage: binauralSH_stop()
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% April 2021

%% ===== Removing Paths ==================================================

% Get the basepath as the directory this function resides in.
% The 'which' solution below is more portable than 'mfilename'
% becase old versions of Matlab does not have "mfilename('fullpath')"
basepath=which('binauralSH_stop');
% Kill the function name from the path.
basepath=basepath(1:end-18);

% Add the base path and the needed sub-directories
f=filesep;
if exist('rmpath')
    rmpath(basepath);
    rmpath([basepath f 'core']);
    rmpath([basepath f 'thirdparty']);
    rmpath([basepath f 'utils']);
end

