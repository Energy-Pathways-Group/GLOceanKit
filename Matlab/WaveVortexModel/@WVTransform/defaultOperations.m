function operations = defaultOperations()
% return array of WVOperation instances initialized by default
%
% This function creates a number of standard StateVariables with associated
% TransformOperations.
%
% - Topic: Internal
% - Declaration: operations = defaultOperations()
% - Returns operations: array of WVOperation instances
operations = WVOperation.empty(0,0);

outputVar(1) = WVVariableAnnotation('Apt',{'k','l','j'},'m/s', 'positive wave coefficients at time t');
outputVar(1).isComplex = 1;

outputVar(2) = WVVariableAnnotation('Amt',{'k','l','j'},'m/s', 'negative wave coefficients at time t');
outputVar(2).isComplex = 1;

outputVar(3) = WVVariableAnnotation('A0t',{'k','l','j'},'m', 'geostrophic coefficients at time t');
outputVar(3).isComplex = 1;

f = @(wvt) wvt.waveVortexCoefficientsAtTimeT();
operations(end+1) = WVOperation('ApAmA0',outputVar,f);

outputVar = WVVariableAnnotation('uMax',{},'m s^{-1}', 'max horizontal fluid speed');
f = @(wvt) max(max(max( sqrt( (wvt.u).^2 + (wvt.v).^2 ) )));
operations(end+1) = WVOperation('uMax',outputVar,f);

outputVar = WVVariableAnnotation('wMax',{},'m s^{-1}', 'max vertical fluid speed');
f = @(wvt) max(max(max( abs(wvt.w)  )));
operations(end+1) = WVOperation('wMax',outputVar,f);

outputVar = WVVariableAnnotation('u',{'x','y','z'},'m/s', 'x-component of the fluid velocity');
outputVar.attributes('standard_name') = 'eastward_sea_water_velocity';
f = @(wvt) wvt.transformToSpatialDomainWithF(wvt.UAp.*wvt.Apt + wvt.UAm.*wvt.Amt + wvt.UA0.*wvt.A0t);
operations(end+1) = WVOperation('u',outputVar,f);

outputVar = WVVariableAnnotation('v',{'x','y','z'},'m/s', 'y-component of the fluid velocity');
outputVar.attributes('standard_name') = 'northward_sea_water_velocity';
f = @(wvt) wvt.transformToSpatialDomainWithF(wvt.VAp.*wvt.Apt + wvt.VAm.*wvt.Amt + wvt.VA0.*wvt.A0t);
operations(end+1) = WVOperation('v',outputVar,f);

outputVar = WVVariableAnnotation('w',{'x','y','z'},'m/s', 'z-component of the fluid velocity');
outputVar.attributes('standard_name') = 'upwardward_sea_water_velocity';
f = @(wvt) wvt.transformToSpatialDomainWithG(wvt.WAp.*wvt.Apt + wvt.WAm.*wvt.Amt);
operations(end+1) = WVOperation('w',outputVar,f);

outputVar = WVVariableAnnotation('p',{'x','y','z'},'kg/m/s2', 'pressure anomaly');
f = @(wvt) wvt.rho0*wvt.g*wvt.transformToSpatialDomainWithF(wvt.NAp.*wvt.Apt + wvt.NAm.*wvt.Amt + wvt.PA0.*wvt.A0t);
operations(end+1) = WVOperation('p',outputVar,f);

outputVar = WVVariableAnnotation('psi',{'x','y','z'},'m^2/s', 'geostrophic streamfunction');
f = @(wvt) wvt.transformToSpatialDomainWithF((wvt.g/wvt.f) * wvt.A0t);
operations(end+1) = WVOperation('psi',outputVar,f);

outputVar = WVVariableAnnotation('eta',{'x','y','z'},'m', 'isopycnal deviation');
f = @(wvt) wvt.transformToSpatialDomainWithG(wvt.NAp.*wvt.Apt + wvt.NAm.*wvt.Amt + wvt.NA0.*wvt.A0t);
operations(end+1) = WVOperation('eta',outputVar,f);

% outputVar = WVVariableAnnotation('rho_e',{'x','y','z'},'kg/m^3', 'excess density');
% f = @(wvt) (wvt.rho0/wvt.g) * shiftdim(wvt.N2,-2) .* wvt.eta;
% operations(end+1) = WVOperation('rho_e',outputVar,f);

outputVar = WVVariableAnnotation('rho_prime',{'x','y','z'},'kg/m3', 'density anomaly');
f = @(wvt) (wvt.rho0/wvt.g) * shiftdim(wvt.N2,-2) .* wvt.eta;
operations(end+1) = WVOperation('rho_prime',outputVar,f);

outputVar = WVVariableAnnotation('rho_total',{'x','y','z'},'kg/m3', 'total potential density');
f = @(wvt) reshape(wvt.rhobar,1,1,[]) + wvt.rho_prime;
operations(end+1) = WVOperation('rho_total',outputVar,f);

outputVar = WVVariableAnnotation('qgpv',{'x','y','z'},'1/s', 'quasigeostrophic potential vorticity');
f = @(wvt) wvt.transformToSpatialDomainWithF( wvt.A0_QGPV_factor .*wvt.A0t);
operations(end+1) = WVOperation('qgpv',outputVar,f);

outputVar = WVVariableAnnotation('seaSurfaceU',{'x','y'},'m/s', 'x-component of the fluid velocity at the surface',detailedDescription='- topic: State Variables');
operations(end+1) = WVOperation('seaSurfaceU', outputVar,@(wvt) wvt.u(:,:,end));

outputVar = WVVariableAnnotation('seaSurfaceV',{'x','y'},'m/s', 'y-component of the fluid velocity at the surface',detailedDescription='- topic: State Variables');
operations(end+1) = WVOperation('seaSurfaceV', outputVar,@(wvt) wvt.v(:,:,end));

outputVar = WVVariableAnnotation('seaSurfaceHeight',{'x','y'},'m', 'sea-surface height');
operations(end+1) = WVOperation('seaSurfaceHeight', outputVar,@(wvt) wvt.p(:,:,end)/(wvt.rho0*wvt.g));

fluxVar(1) = WVVariableAnnotation('Fp',{'k','l','j'},'m/s2', 'non-linear flux into Ap', isComplex=1, detailedDescription='- topic: State Variables');
fluxVar(2) = WVVariableAnnotation('Fm',{'k','l','j'},'m/s2', 'non-linear flux into Am', isComplex=1,detailedDescription='- topic: State Variables');
fluxVar(3) = WVVariableAnnotation('F0',{'k','l','j'},'m/s', 'non-linear flux into A0', isComplex=1,detailedDescription='- topic: State Variables');
operations(end+1) = WVOperation('nonlinearFlux',fluxVar,@(wvt) wvt.nonlinearFlux);

end