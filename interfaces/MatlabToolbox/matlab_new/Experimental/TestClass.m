classdef TestClass < matlab.unittest.TestCase
    
    properties (Constant = true)
        tolerance = 1e-9;
    end
    
    %% Test Method Block
    
    methods (Test)
        
        %% Test Functions
        
        function testClassConstructor(testCase)
            % Verify that the class constructor creates all the properties
            % necessary for ThermoPhase class
            
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            allProperties = properties(tpTest)';
            testCase.verifyEqual(allProperties, PropList);
        end
        
        function testPropertyT(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{2});
            ref = ThermoRef{2};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyP(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{3});
            ref = ThermoRef{3};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyD(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{4});
            ref = ThermoRef{4};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyX(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{5});
            ref = ThermoRef{5};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyY(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{6});
            ref = ThermoRef{6};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyH(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{7});
            ref = ThermoRef{7};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyS(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{8});
            ref = ThermoRef{8};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyU(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{9});
            ref = ThermoRef{9};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyG(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{10});
            ref = ThermoRef{10};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyV(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{11});
            ref = ThermoRef{11};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertybasis(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{12});
            ref = ThermoRef{12};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyDP(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{13});
            ref = ThermoRef{13};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyDPX(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{14});
            ref = ThermoRef{14};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyDPY(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{15});
            ref = ThermoRef{15};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyHP(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{16});
            ref = ThermoRef{16};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyHPX(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{17});
            ref = ThermoRef{17};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyHPY(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{18});
            ref = ThermoRef{18};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyPV(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{19});
            ref = ThermoRef{19};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyPVX(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{20});
            ref = ThermoRef{20};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyPVY(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{21});
            ref = ThermoRef{21};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertySH(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{22});
            ref = ThermoRef{22};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertySHX(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{23});
            ref = ThermoRef{23};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertySHY(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{24});
            ref = ThermoRef{24};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertySP(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{25});
            ref = ThermoRef{25};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertySPX(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{26});
            ref = ThermoRef{26};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertySPY(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{27});
            ref = ThermoRef{27};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyST(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{28});
            ref = ThermoRef{28};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertySTX(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{29});
            ref = ThermoRef{29};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertySTY(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{30});
            ref = ThermoRef{30};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertySV(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{31});
            ref = ThermoRef{31};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertySVX(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{32});
            ref = ThermoRef{32};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertySVY(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{33});
            ref = ThermoRef{33};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyTD(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{34});
            ref = ThermoRef{34};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyTDX(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{35});
            ref = ThermoRef{35};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyTDY(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{36});
            ref = ThermoRef{36};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyTH(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{37});
            ref = ThermoRef{37};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyTHX(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{38});
            ref = ThermoRef{38};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyTHY(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{39});
            ref = ThermoRef{39};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyTP(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{40});
            ref = ThermoRef{40};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyTPX(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{41});
            ref = ThermoRef{41};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyTPY(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{42});
            ref = ThermoRef{42};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyTV(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{43});
            ref = ThermoRef{43};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyTVX(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{44});
            ref = ThermoRef{44};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyTVY(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{45});
            ref = ThermoRef{45};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyUV(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{46});
            ref = ThermoRef{46};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyUVX(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{47});
            ref = ThermoRef{47};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyUVY(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{48});
            ref = ThermoRef{48};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyUP(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{49});
            ref = ThermoRef{49};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyUPX(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{50});
            ref = ThermoRef{50};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyUPY(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{51});
            ref = ThermoRef{51};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyVH(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{52});
            ref = ThermoRef{52};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyVHX(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{53});
            ref = ThermoRef{53};
            testCase.verifyEqual(prop, ref);
        end

        function testPropertyVHY(testCase)
            tpTest = ThermoPhase('gri30.yaml', 'gri30');
            load('ThermoPhaseTest.mat');
            prop = tpTest.(PropList{54});
            ref = ThermoRef{54};
            testCase.verifyEqual(prop, ref);
        end
       
    end
    
end