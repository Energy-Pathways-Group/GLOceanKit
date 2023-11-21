%% Initialize a wave vortex transform

N0 = 3*2*pi/3600;
L_gm = 1300;
N2 = @(z) N0*N0*exp(2*z/L_gm);
wvt = WVTransformHydrostatic([800e3, 400e3, 4000],[128, 64, 40], N2=N2,latitude=30);
% omega = wvt.initWithWaveModes(k=10,l=0,j=1,phi=0,U=0.2,sign=1);
% period = 2*pi/omega;
wvt.initWithGMSpectrum(1.0);
period = wvt.inertialPeriod;
%% Initialize a WVModel (model) with the WVTransform
% You can also (optionally) specify a nonlinear flux operation. Here we add 
% the nonlinear Boussinesq flux.

model = WVModel(wvt,nonlinearFlux=WVNonlinearFlux(wvt,shouldAntialias=1,uv_damp=wvt.uMax));
%% Add particles with float like behavior
% _Floats_ will faithfully follow the fluid flow, while _drifters_ will only 
% follow (u,v), staying a fixed depth.

nTrajectories = 101;
xFloat = wvt.Lx/2*ones(1,nTrajectories);
yFloat = wvt.Ly/2*ones(1,nTrajectories);
zFloat = linspace(-wvt.Lz,0,nTrajectories);
model.setFloatPositions(xFloat,yFloat,zFloat);
%% Setup the model integrator, integrate 3 wave periods
% Tell the model the size of the time step we want to take as well as output 
% interval. Often the integrator can run using only an advective time-step constraint, 
% but when advecting particles, it is important to use the oscillatory time step, 
% which is typically much shorter.
% 
% -setupIntegrator will automatically be called if you do not call it.

model.setupIntegrator(timeStepConstraint="oscillatory", outputInterval=period/10);
model.createNetCDFFileForModelOutput('PlaneWaveWithFloats.nc',shouldOverwriteExisting=1);
model.integrateToTime(3*wvt.inertialPeriod);
%% Read the particles from file

ncfile = model.ncfile;
[x,y,z] = ncfile.readVariables('float_x','float_y','float_z');

figure, plot3(x.',y.',z.')