% Use these parameters to generate a few synthetic trajectories.
sigma = 8e-6;
theta = 30*pi/180;
zeta = 0;
kappa = 1;

velocityField = LinearVelocityField(sigma,theta,zeta);
integrator = AdvectionDiffusionIntegrator(velocityField,kappa);

% Particles placed in a "plus" configuration.
x0 = [-500; -250; 0; 0; 0; 0; 0; 250; 500];
y0 = [0; 0; -500; -250; 0; 250; 500; 0; 0;];

% Now compute the trajectories
T = 3*86400;
dt = 4*3600;
[t,x,y] = integrator.particleTrajectories(x0,y0,T,dt);

% And then try to estimate the parameters 
parametersToEstimate = [ModelParameter.strain];
parameterEstimates = EstimateLinearVelocityFieldParameters( x, y, t, parametersToEstimate );