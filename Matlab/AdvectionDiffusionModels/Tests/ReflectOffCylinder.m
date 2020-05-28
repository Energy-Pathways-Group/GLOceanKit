%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Cylinder flow
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

model = CylinderFlow();
% figure
% model.plotStreamfunction(), hold on
% model.plotVelocityField();

kappa = 1e3;
integrator = AdvectionDiffusionIntegrator(model,kappa);

% determine reasonable integration time scales
dt = 864;

sigma_dx = sqrt(2*kappa*dt);
fprintf('Typical particle step is %.2f km\n',sigma_dx*1e-3)

% place particles throughout the valid domain
x = linspace(min(model.xlim),max(model.xlim),10);
y = linspace(min(model.ylim),max(model.ylim),10);
[x0,y0] = ndgrid(x,y);

[t,x,y] = integrator.particleTrajectories(x0,y0,dt*100,dt);

figure
model.plotVelocityField(), hold on
model.plotTrajectories(x,y)
