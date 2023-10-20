function variableAnnotations = defaultVariableAnnotations()
% return array of WVVariableAnnotation instances initialized by default
%
% This function creates annotations for the built-in variables supported by
% the WVTransform.
%
% - Topic: Internal
% - Declaration: operations = defaultVariableAnnotations()
% - Returns operations: array of WVVariableAnnotation instances
variableAnnotations = WVVariableAnnotation.empty(0,0);

variableAnnotations(end+1) = WVVariableAnnotation('t',{}, 's', 'time of observations');

annotation = WVVariableAnnotation('A0',{'k','l','j'},'m', 'geostrophic coefficients at reference time t0');
annotation.isComplex = 1;
annotation.isVariableWithLinearTimeStep = 0;
annotation.isVariableWithNonlinearTimeStep = 1;
variableAnnotations(end+1) = annotation;

annotation = WVVariableAnnotation('Ap',{'k','l','j'},'m/s', 'positive wave coefficients at reference time t0');
annotation.isComplex = 1;
annotation.isVariableWithLinearTimeStep = 0;
annotation.isVariableWithNonlinearTimeStep = 1;
variableAnnotations(end+1) = annotation;

annotation = WVVariableAnnotation('Am',{'k','l','j'},'m/s', 'negative wave coefficients at reference time t0');
annotation.isComplex = 1;
annotation.isVariableWithLinearTimeStep = 0;
annotation.isVariableWithNonlinearTimeStep = 1;
variableAnnotations(end+1) = annotation;

annotation = WVVariableAnnotation('totalEnergy',{},'m3/s2', 'horizontally-averaged depth-integrated energy computed spectrally from wave-vortex coefficients');
annotation.isVariableWithLinearTimeStep = 0;
annotation.isVariableWithNonlinearTimeStep = 1;
variableAnnotations(end+1) = annotation;

annotation = WVVariableAnnotation('totalEnergySpatiallyIntegrated',{},'m3/s2', 'horizontally-averaged depth-integrated energy computed in the spatial domain');
annotation.isVariableWithLinearTimeStep = 0;
annotation.isVariableWithNonlinearTimeStep = 1;
variableAnnotations(end+1) = annotation;

annotation = WVVariableAnnotation('totalHydrostaticEnergy',{},'m3/s2', 'horizontally-averaged depth-integrated energy *without w* computed in the spatial domain');
annotation.isVariableWithLinearTimeStep = 0;
annotation.isVariableWithNonlinearTimeStep = 1;
variableAnnotations(end+1) = annotation;

annotation = WVVariableAnnotation('internalWaveEnergyPlus',{},'m3/s2', 'total energy, internal waves, positive');
annotation.isVariableWithLinearTimeStep = 0;
annotation.isVariableWithNonlinearTimeStep = 1;
variableAnnotations(end+1) = annotation;

annotation = WVVariableAnnotation('internalWaveEnergyMinus',{},'m3/s2', 'total energy, internal waves, minus');
annotation.isVariableWithLinearTimeStep = 0;
annotation.isVariableWithNonlinearTimeStep = 1;
variableAnnotations(end+1) = annotation;

annotation = WVVariableAnnotation('waveEnergy',{},'m3/s2', 'total energy, waves');
annotation.isVariableWithLinearTimeStep = 0;
annotation.isVariableWithNonlinearTimeStep = 1;
variableAnnotations(end+1) = annotation;

annotation = WVVariableAnnotation('inertialEnergyBaroclinic',{},'m3/s2', 'total energy, inertial oscillations, baroclinic');
annotation.isVariableWithLinearTimeStep = 0;
annotation.isVariableWithNonlinearTimeStep = 1;
variableAnnotations(end+1) = annotation;

annotation = WVVariableAnnotation('inertialEnergyBarotropic',{},'m3/s2', 'total energy, inertial oscillations, barotropic');
annotation.isVariableWithLinearTimeStep = 0;
annotation.isVariableWithNonlinearTimeStep = 1;
variableAnnotations(end+1) = annotation;

annotation = WVVariableAnnotation('inertialEnergy',{},'m3/s2', 'total energy, inertial oscillations');
annotation.isVariableWithLinearTimeStep = 0;
annotation.isVariableWithNonlinearTimeStep = 1;
variableAnnotations(end+1) = annotation;

annotation = WVVariableAnnotation('geostrophicEnergyBaroclinic',{},'m3/s2', 'total energy, geostrophic, baroclinic');
annotation.isVariableWithLinearTimeStep = 0;
annotation.isVariableWithNonlinearTimeStep = 1;
variableAnnotations(end+1) = annotation;

annotation = WVVariableAnnotation('geostrophicEnergyBarotropic',{},'m3/s2', 'total energy, geostrophic, barotropic');
annotation.isVariableWithLinearTimeStep = 0;
annotation.isVariableWithNonlinearTimeStep = 1;
variableAnnotations(end+1) = annotation;

annotation = WVVariableAnnotation('geostrophicEnergy',{},'m3/s2', 'total energy, geostrophic');
annotation.isVariableWithLinearTimeStep = 0;
annotation.isVariableWithNonlinearTimeStep = 1;
variableAnnotations(end+1) = annotation;

end