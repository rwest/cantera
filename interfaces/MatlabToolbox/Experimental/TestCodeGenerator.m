load('ThermoPhaseTest.mat');

fid = fopen('Experimental/code.m', 'w');

for i = 2:length(PropList)
    generatePropertyGetTestMethod(fid, i, PropList);
end

fclose(fid);

clear all
close all

function generatePropertyGetTestMethod(fileID, i, propList)
    % Generate code for testing properties
    fprintf(fileID, 'function testProperty%s%s\r', propList{i}, ...
            '(testCase)');
    fprintf(fileID, '    tpTest = testCase.tp;\r');
    fprintf(fileID, '    prop = tpTest.(testCase.propList{%d});\r', i);
    fprintf(fileID, '    ref = testCase.refDB{%d};\r', i);
    fprintf(fileID, '    testCase.verifyEqual(prop, ref);\r');
    fprintf(fileID, 'end\r');
    fprintf(fileID, '\r');
end