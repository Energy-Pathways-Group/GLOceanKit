%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Divergence/Convergence box
%
% Here we start with a bunch of particles evenly spread throughout the
% domain, but in this case, there's a weak convergence in half the box and
% a weak divergence in the other half of the box.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

showEndStateOnly = 1;
FramesFolder = './FramesScratch';
if exist(FramesFolder,'dir') == 0
	mkdir(FramesFolder);
end

model = DivergenceBox();

kappa = 100;
integrator = AdvectionDiffusionIntegrator(model,kappa);

% determine reasonable integration time scales
T = 2*86400;
dt = 864;

% place particles on the left-quarter of the domain
x = linspace(min(model.xlim),max(model.xlim),20);
y = linspace(min(model.ylim),max(model.ylim),10);
[x0,y0] = ndgrid(x,y);

n = floor(T/dt/3);
[t1,x1,y1] = integrator.particleTrajectories(x0,y0,n*dt,dt);
[t2,x2,y2] = integrator.particleTrajectories(x1(end),y1(end),n*dt,dt);
[t3,x3,y3] = integrator.particleTrajectories(x2(end),y2(end),n*dt,dt);


figure('Position', [50 50 680 400])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');


if showEndStateOnly == 1
    iTime0 = length(t);
else
    iTime0 = 1;
end
for iTime=iTime0:length(t)
    clf
    
    model.plotBounds(), hold on
    divbox.plotVelocityField
    axis equal
    xticks([])
    yticks([])
    title([])
    axis off
    
    scatter( model.visualScale*x(iTime,:), model.visualScale*y(iTime,:), 8^2, 'r', 'fill')
    
    text(0,model.visualScale*model.Ly+.2,'Weak divergence')
    
    if showEndStateOnly == 1
        continue;
    end
    
    % write everything out
    output = sprintf('%s/t_%03d', FramesFolder,iTime-1);
    print('-depsc2', output)
end