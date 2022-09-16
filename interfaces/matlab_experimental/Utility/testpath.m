% Set up environment for testing the Cantera Matlab interface
% from within the Cantera source tree. Run this file from the
% root of the Cantera source tree, for example:
%
%    cd ~/src/cantera
%    run interfaces/matlab/testpath.m

% Unload the mex file so copying the DLL will work
clear ctmethods

% get the list of directories on the Matlab path
dirs = regexp(path, ['([^' pathsep ']*)'], 'match');

% if 'cantera' is already in the path, remove it
for i = 1:length(dirs)
    if strfind(dirs{i}, 'Cantera')
        rmpath(dirs{i});
        continue;
    end
    if strfind(dirs{i}, 'cantera')
        rmpath(dirs{i});
    end
end

% Add the Cantera toolbox to the Matlab path
path(path, genpath(fullfile(pwd, '..')));
cantera_root = fullfile(pwd, '..', '..', '..');

% Set path to Python module
if strcmp(getenv('PYTHONPATH'), '')
    setenv('PYTHONPATH', fullfile(cantera_root, 'build', 'python'))
end

% A simple test to make sure that the ctmethods.mex file is present and working
LoadCantera;

f = Func('polynomial', 3, [1,2,3,4]);
if f(1) == 10
    disp('Cantera is successfully loaded.')
else
    disp('Something is wrong with the LoadCantera.')
end

adddir(fullfile(cantera_root, 'data', 'inputs'))
