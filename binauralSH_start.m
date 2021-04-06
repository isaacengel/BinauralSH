function binauralSH_start(showinfo)
%binauralSH_start Start the Binaural SH Toolbox
%
%   Usage: binauralSH_start(showinfo)
%
% INPUT:
%   showinfo = 0: print nothing (default)
%              1: print version and web page link
%
% AUTHOR: Isaac Engel - isaac.engel(at)imperial.ac.uk
% April 2021

%% ===== Checking of input  parameters ===================================
if ~exist('showinfo','var')
    showinfo = 0;
end

%% ===== Adding Paths ====================================================

% Get the basepath as the directory this function resides in.
% The 'which' solution below is more portable than 'mfilename'
% becase old versions of Matlab does not have "mfilename('fullpath')"
basepath=which('binauralSH_start');
% Kill the function name from the path.
basepath=basepath(1:end-19);

% Add the base path and the needed sub-directories
f=filesep;
if exist('addpath')
    addpath(basepath);
    addpath([basepath f 'core']);
    addpath([basepath f 'examples']);
    addpath([basepath f 'hrtfs']);
    addpath(genpath([basepath f 'thirdparty']));
    addpath([basepath f 'utils']);
else
    path(path,basepath);
    path(path,[basepath f 'core']);
    path(path,[basepath f 'examples']);
    path(path,[basepath f 'hrtfs']);
    path(path,genpath([basepath f 'thirdparty']));
    path(path,[basepath f 'utils']);
end

%% ===== Banner ==========================================================
if showinfo
    disp(['Binaural SH tooolbox version ' binauralSH_version '.']);
    disp('Github: https://github.com/isaacengel/BinauralSH')
    disp('Contact: Isaac Engel - isaac.engel(at)imperial.ac.uk') 
end

