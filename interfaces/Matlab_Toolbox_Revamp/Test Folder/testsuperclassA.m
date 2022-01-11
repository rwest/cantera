classdef testsuperclassA < handle
    properties 
        T
        P
        tp_id
    end
    methods
        function tp = testsuperclassA(temperature, pressure)
            tp.T = temperature;
            tp.P = pressure;
            tp.tp_id = 1;
        end
        % get methods
        function temperature = get.T(tp)
            temperature = tp.T;
        end
        function pressure = get.P(tp)
            pressure = tp.P;
        end
        % set methods
        function tp = set.T(tp, temperature)
            tp.T = temperature;
        end
        function tp = set.P(tp, pressure)
            tp.P = pressure;
        end
    end
end