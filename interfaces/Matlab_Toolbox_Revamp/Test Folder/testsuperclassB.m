classdef testsuperclassB < handle
    properties
        mdot
        rho
        tr_id
    end
    methods
        function tr = testsuperclassB(flowrate, density)
            tr.mdot = flowrate;
            tr.rho = density;
            tr.tr_id = 2;
        end
    end
end