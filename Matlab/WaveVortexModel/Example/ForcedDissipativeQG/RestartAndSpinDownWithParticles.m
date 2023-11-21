%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize a new WVModel (which includes a WVTransform
% and a WVNonlinearFluxOperation) from existing output. We will start from
% the final time-point.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Loading the model from file assumes that you want to use the same forcing
% as before...
model = WVModel.modelFromFile('ForcedDissipativeQG-spinup-256.nc',restartIndex=Inf);

%...so let's replace that forcing with an unforced version, but still use
%the same damping (nu).
unforcedFlux = QGPVE(model.wvt,nu=model.nonlinearFluxOperation.nu);
model.nonlinearFluxOperation = unforcedFlux;

[xFloat,yFloat] = ndgrid(wvt.x(1:4:end),wvt.y(1:4:end));
xFloat = reshape(xFloat,1,[]);
yFloat = reshape(yFloat,1,[]);
nTrajectories = length(xFloat);
model.setDrifterPositions(xFloat,yFloat,[],'zeta','eta_f');


% Setting the output interval based on the enstrophy time scale
model.setupIntegrator(deltaT=0.5*unforcedFlux.dampingTimeScale,outputInterval=30*60);
model.createNetCDFFileForModelOutput('ForcedDissipativeSpinDownQG-particles-256.nc',shouldOverwriteExisting=1);
model.setNetCDFOutputVariables('A0','psi','zeta_z');

model.integrateToTime(wvt.t + 5*86400);

ncfile = model.ncfile;
[x,y] = ncfile.readVariables('drifter_x','drifter_y');
figure, plot(x.',y.')