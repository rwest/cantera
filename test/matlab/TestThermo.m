classdef TestThermo < TestCase
    properties
        gas
    end

    methods
        function self = TestThermo(name)
            self = self@TestCase(name);
        end

        function setUp(self)
            global staticTestThermoGas
            if isempty(staticTestThermoGas)
                staticTestThermoGas = Solution('../data/steam-reforming.yaml', 'full');
            end
            self.gas = staticTestThermoGas;
            gas.TPY = {300, oneatm, [0.5, 0, 0.5, 0, 0, 0, 0]};
        end

%        function tearDown(self)
%        end

        function testCounts(self)
            assertEqual(self.gas.nElements, 4)
            assertEqual(self.gas.nSpecies, 7)
        end

        function testElements(self)
            for i = 1:nElements(self.gas)
                name = self.gas.elementName(i);
                assertEqual(i, self.gas.elementIndex(name))
            end
        end

        function testSpecies(self)
            for i = 1:self.gas.nSpecies
                name = self.gas.speciesName(i);
                assertEqual(i, self.gas.speciesIndex(name))
            end
        end

        function test_nAtoms(self)
           assertEqual(self.gas.nAtoms(1, 3), 4)
           assertEqual(self.gas.nAtoms(1, 4), 0)
           assertEqual(self.gas.nAtoms(3, 2), 2)
           assertExceptionThrown(@() self.gas.nAtoms(2, 5), '')
           assertExceptionThrown(@() self.gas.nAtoms(8, 1), '')
        end

        function testSetState(self)
            self.gas.basis = 'mass';
            u0 = self.gas.U;
            h0 = self.gas.H;
            s0 = self.gas.S;
            v0 = 1/self.gas.D;
            T0 = self.gas.T;
            P0 = self.gas.P;

            set(self.gas, 'T', 400, 'P', 5*oneatm);
            gas.TP = {400, 5*oneatm};
            assertAlmostEqual(self.gas.T, 400)
            assertAlmostEqual(self.gas.P, 5*oneatm)

            gas.HP = {h0, P0};
            assertAlmostEqual(self.gas.T, T0, 1e-8)
            assertAlmostEqual(self.gas.S, s0, 1e-8)

            gas.TP = {400, 5*oneatm};
            gas.UV = {u0, v0};
            assertAlmostEqual(self.gas.P, P0, 1e-8)
            assertAlmostEqual(self.gas.H, h0, 1e-8)
        end
    end
end
