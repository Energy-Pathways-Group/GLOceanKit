%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Simple box
%
% The idea here is that we start a bunch of particles in part of the box
% and after some time, the particles evenly spread to all parts of the box.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

showEndStateOnly = 0;
FramesFolder = './FramesScratch';
if exist(FramesFolder,'dir') == 0
	mkdir(FramesFolder);
end

box = SimpleBox();

kappa = 100;
integrator = AdvectionDiffusionIntegrator(box,kappa);

% determine reasonable integration time scales
T = 2*86400;
dt = 864;

% place particles on the left-quarter of the domain
x = linspace(min(box.xlim),max(box.xlim)/4,20);
y = linspace(min(box.ylim),max(box.ylim),10);
[x0,y0] = ndgrid(x,y);

[t,x,y] = integrator.particleTrajectories(x0,y0,T,dt);

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
    
    box.plotBounds(), hold on
    axis equal
    xticks([])
    yticks([])
    axis off
    
    scatter( box.visualScale*x(iTime,:), box.visualScale*y(iTime,:), 8^2, 'r', 'fill')
    
    if showEndStateOnly == 1
        continue;
    end
    
    % write everything out
    output = sprintf('%s/t_%03d', FramesFolder,iTime-1);
    print('-depsc2', output)
end