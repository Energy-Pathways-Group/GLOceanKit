%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Meandering Jet
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

jet = MeanderingJet();

kappa = 0;
integrator = AdvectionDiffusionIntegrator(jet,kappa);

% determine reasonable integration time scales
T = 100*jet.Lx/jet.U;
dt = 864;

% place particles in a row along the jet
x0 = zeros(100,1);
y0 = linspace(0,100e3,100);

[t,x,y] = integrator.particleTrajectories(x0,y0,T,dt);

figure
jet.plotVelocityField(), hold on
jet.plotTrajectories(x,y)
print('figures/jet_trajectories_kappa0.png','-dpng')

cv = diff(x+sqrt(-1)*y)/dt;
[f,spp,snn,spn] = mspec(dt,cv,[]);
figure

subplot(2,1,1)
plot(f,spp), hold on
plot(f,median(spp,2),'k','LineWidth',2)
vlines( 2*pi/(jet.Lx/jet.U) )
xlim([min(f) 1e-3])
ylim([1e0 1e7])
xlog, ylog
subplot(2,1,2)
plot(f,snn), hold on
plot(f,median(snn,2),'k','LineWidth',2)
vlines( 2*pi/(jet.Lx/jet.U) )
xlog, ylog
xlim([min(f) 1e-3])
ylim([1e0 1e7])
print('figures/jet_spectra_kappa0.png','-dpng')