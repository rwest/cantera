clear all
close all
cleanup

a = ThermoPhase('gri30.yaml', 'gri30');
PropList = properties(a)';
MethodList = methods(a)';

for i = 2:length(PropList)
    ThermoRef{i} = a.(PropList{i});
end

save('Experimental/ThermoPhaseTest.mat', 'PropList', 'ThermoRef');
