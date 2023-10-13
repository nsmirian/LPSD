function mydir = my_startup()

restoredefaultpath;
rehash toolboxcache

ones(10)*ones(10);
fspecial('gaussian');

me = mfilename;
mydir = which(me);
mydir = mydir(1:end-2-numel(me));
addpath(mydir);
% addpath(genpath([mydir 'Tools']));
% addpath(genpath([mydir 'Utils']));
% addpath(genpath([mydir 'calib_vals']));
addpath('/home/ttflinac/bin/matlab')
addpath('/local/lib')
addpath('/home/ttflinac/released_software/matlab/hlc_toolbox_common/')

addpath('/home/kammer/MATLAB/createElogEntry/');