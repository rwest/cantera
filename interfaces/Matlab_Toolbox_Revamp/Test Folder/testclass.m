classdef testclass < handle
    properties 
        T
        P
    end
    methods
        function tp = testclass(temperature, pressure)
            tp.T = temperature;
            tp.P = pressure;
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