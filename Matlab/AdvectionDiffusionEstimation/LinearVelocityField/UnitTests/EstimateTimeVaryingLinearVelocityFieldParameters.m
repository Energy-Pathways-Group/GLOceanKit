%
% Meandering Jet
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% sigma = 4e-6;
% kappa = 0.5;
% T = round(1/sigma/86400)*86400;
% T = 15*86400;
% dt = 2*3600;
% 
% velocityField = LinearVelocityField(sigma,0,0);
% integrator = AdvectionDiffusionIntegrator(velocityField,kappa);
% 
% x0 = [-500; -250; 0; 0; 0; 0; 0; 250; 500];
% y0 = [0; 0; -500; -250; 0; 250; 500; 0; 0;];
% 
% 
% [t,x,y] = integrator.particleTrajectories(x0,y0,T,dt);
% dof = 2;
% K = 2;
[parameters,error,B] = FitTrajectoriesToTimeVaryingEllipseModel( x, y, t, 'strain-diffusive', dof, K );

figure
subplot(1,3,1)
plot(t/86400,B*parameters.m_sigma)
subplot(1,3,2)
plot(t/86400,B*parameters.m_theta*180/pi)
subplot(1,3,3)
plot(t/86400,B*parameters.m_kappa)