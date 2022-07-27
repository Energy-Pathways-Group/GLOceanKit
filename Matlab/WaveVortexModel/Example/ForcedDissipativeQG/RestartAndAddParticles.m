%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize a new WaveVortexModel (which includes a WaveVortexTransform
% and a WVNonlinearFluxOperation) from existing output. We will start from
% the final time-point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = WaveVortexModel.modelFromFile('ForcedDissipativeQG-spinup-512.nc',restartIndex=Inf);

wvt = model.wvt;
outputVar = WVVariableAnnotation('eta',{'x','y','z'},'1/s^3', 'enstrophy forcing');
f = @(wvt) wvt.zeta_z .* wvt.F_psi;
wvt.addOperation(WVOperation('eta',outputVar,f));

[xFloat,yFloat] = ndgrid(wvt.x(1:4:end),wvt.y(1:4:end));
xFloat = reshape(xFloat,1,[]);
yFloat = reshape(yFloat,1,[]);
nTrajectories = length(xFloat);
model.setDrifterPositions(xFloat,yFloat,[],'zeta_z','eta');


% Setting the output interval based on the enstrophy time scale
model.setupIntegrator(deltaT=0.5*model.nonlinearFlux.dampingTimeScale,outputInterval=30*60);
model.createNetCDFFileForModelOutput('ForcedDissipativeQG-particles-512.nc',shouldOverwriteExisting=1);
model.setNetCDFOutputVariables('A0','psi','zeta_z','F_psi','F0_psi');

model.integrateToTime(wvt.t + 10*86400);

ncfile = model.ncfile;
[x,y] = ncfile.readVariables('drifter-x','drifter-y');
figure, plot(x.',y.')