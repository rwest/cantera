function  runCanteraTests(varargin)
% Set Cantera paths
run ../../interfaces/matlab/Utility/testpath.m

% Set test framework path
path('../../ext/matlab_xunit', path)

% run the tests
if nargin == 0
    runtests -verbose
else
    runtests('-verbose','-logfile',varargin{1})
end

% unload the Cantera Clib
cleanup
UnloadCantera

% delete global objects created by some of the test
clear global static*