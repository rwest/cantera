classdef docuTest 
    properties
        propA
        propB
        id
    end
    
    methods
        function d = docuTest(a)
            % DOCUTEST docuTest class constructor.
            % d = docuTest(a)
            % :param a: Arbitrary number.
            % :return:
            %    Instance of class :mat:func:`docuTest`
            d.id = 1;
            d.propA = 42;
            d.propB = 13;
        end
        
        function v = aplusb(d)
            % APLUSB Return sum of property A and B
            % v = d.aplusb
            % :param d:
            %    Instance of class :mat:func:`docuTest`
            v = d.propA + d.propB;
        end
    end
        
end