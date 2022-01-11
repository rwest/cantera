classdef testclassA < handle & testsuperclassA & testsuperclassB
    methods
        function t = testclassA(temperature, pressure, flowrate, density)
            t@testsuperclassA(temperature, pressure);
            t@testsuperclassB(flowrate, density);
        end
    end
end