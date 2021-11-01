function testfunc(varargin)
    disp("number of inputs: " + nargin)
    celldisp(varargin(2:end))
end