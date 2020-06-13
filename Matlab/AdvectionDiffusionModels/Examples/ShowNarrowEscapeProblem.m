model = NarrowEscapeProblem();

kappa = 100;
integrator = AdvectionDiffusionIntegrator(model,kappa);

% determine reasonable integration time scales
T = 2*86400;
dt = 864;

% place particles inside the box
x = linspace(0,model.L,20);
y = linspace(0,model.L,20);
[x0,y0] = ndgrid(x,y);

[t,x,y] = integrator.particleTrajectories(x0,y0,T,dt);

figure
plot(scale(model.obstacles,1/model.visualScale),'FaceColor','black','FaceAlpha',1.0)
axis equal
hold on
scatter(x0(:)/model.visualScale,y0(:)/model.visualScale,'filled')
plot(x/model.visualScale,y/model.visualScale)