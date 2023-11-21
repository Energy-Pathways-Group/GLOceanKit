%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize a new WVModel (which includes a WVTransform
% and a WVNonlinearFluxOperation) from existing output. We will start from
% the final time-point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = WVModel.modelFromFile('ForcedDissipativeQG-spinup-512.nc',restartIndex=Inf);

wvt = model.wvt;
outputVar = WVVariableAnnotation('eta_f',{'x','y','z'},'1/s^3', 'enstrophy forcing');
f = @(wvt) wvt.qgpv .* wvt.F_psi;
wvt.addOperation(WVOperation('eta_f',outputVar,f));

[xFloat,yFloat] = ndgrid(wvt.x(1:4:end),wvt.y(1:4:end));
xFloat = reshape(xFloat,1,[]);
yFloat = reshape(yFloat,1,[]);
nTrajectories = length(xFloat);
model.setDrifterPositions(xFloat,yFloat,[],'qgpv','eta');


% Setting the output interval based on the enstrophy time scale
model.setupIntegrator(deltaT=0.5*model.nonlinearFluxOperation.dampingTimeScale,outputInterval=30*60);
model.createNetCDFFileForModelOutput('ForcedDissipativeQG-particles-512.nc',shouldOverwriteExisting=1);
model.setNetCDFOutputVariables('A0','psi','qgpv','F_psi','F0_psi');

model.integrateToTime(wvt.t + 10*86400);

ncfile = model.ncfile;
[x,y] = ncfile.readVariables('drifter_x','drifter_y');
figure, plot(x.',y.')