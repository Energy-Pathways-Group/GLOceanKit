%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% NarrowEscapeProblem with Star
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

showEndStateOnly = 0;
FramesFolder = './FramesScratch';
if exist(FramesFolder,'dir') == 0
	mkdir(FramesFolder);
end

model = NarrowEscapeProblem();

R = 1e3;
% theta = linspace(0,2*pi-2*pi/30,30);
% model.obstacles = cat(1,model.obstacles, polyshape(R*cos(theta)+model.L/2,R*sin(theta)+model.L/2));

np = 5;
x = zeros(np,1);
y = zeros(np,1);
theta = linspace(0,2*pi-2*pi/(2*np),(2*np));
for i=1:(2*np)
    if mod(i,2) == 1
        x(i) = 2*R*cos(theta(i))+model.L/2;
        y(i) = 2*R*sin(theta(i))+model.L/2;
    else
        x(i) = 0.6*R*cos(theta(i))+model.L/2;
        y(i) = 0.6*R*sin(theta(i))+model.L/2;
    end
end
model.obstacles = cat(1,model.obstacles, polyshape(x,y));

kappa = 20;
integrator = AdvectionDiffusionIntegrator(model,kappa);

maxMovieLength = 20; % seconds
totalMovieFrames = maxMovieLength*30;

% determine reasonable integration time scales
T = 2*86400;
dt = T/totalMovieFrames;

% place particles inside the box
x = linspace(0,model.L,20);
y = linspace(0,model.L,20);
[x0,y0] = ndgrid(x,y);

[t,x,y] = integrator.particleTrajectories(x0,y0,T,dt);

fig = figure('Position', [50 50 800 800]);
set(gcf,'PaperPositionMode','manual')
set(gcf, 'Color', 'w');


if showEndStateOnly == 1
    iTime0 = length(t);
else
    iTime0 = 1;
end

tailLength=15;
for iTime=iTime0:length(t)
    clf
    
    plot(scale(model.obstacles,1/model.visualScale),'FaceColor','black','FaceAlpha',1.0); hold on
    
    startIndex = iTime-tailLength;
    if (startIndex < 1) startIndex=1; end
    range=startIndex:iTime;
    
    model.plotTrajectories(x(range,:),y(range,:),'LineWidth',1.5)
    
    text(0,(model.L+1.5*model.delta)/model.visualScale,sprintf('total escapees: %d', sum(x(range(end),:)<0)));
    
    axis equal
    xlim(model.xVisLim/model.visualScale);
    ylim(model.xVisLim/model.visualScale);
    xticks([])
    yticks([])
    axis off
    
    plt = gca;
    plt.Position = [0 0 1 1];
    
    if showEndStateOnly == 1
        continue;
    end
    
    % write everything out
    output = sprintf('%s/t_%03d', FramesFolder,iTime-1);
    print('-dpng','-r300', output)
end

% figure
% plot(scale(model.obstacles,1/model.visualScale),'FaceColor','black','FaceAlpha',1.0)
% axis equal
% hold on
% scatter(x0(1,:)/model.visualScale,y0(1,:)/model.visualScale,'filled')
% plot(x/model.visualScale,y/model.visualScale)