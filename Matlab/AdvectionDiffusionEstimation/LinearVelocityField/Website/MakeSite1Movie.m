% Load our giant list of possible paths.
SiteNumber=1;
K = 4;
dof = 4;
iModel = 3;
totalPermutations = 1000;
shouldSaveFigures = 1;
framesFolder = 'Site1Frames';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Make the frames folder
%
if exist(framesFolder, 'dir') == 0
	mkdir(framesFolder);
end

load(sprintf('../FluidsPaperFigures/smoothedGriddedRho%dDrifters.mat',SiteNumber));
% The last drifter at both sites is only partial time series.
x = x(:,1:(end-1));
y = y(:,1:(end-1));

filename = sprintf('../FluidsPaperFigures/BootstrapData/Rho%dDrifterSplineFits1000_dof%d.mat',SiteNumber,dof);
load(filename);
[~,mostLikelyIndices] = sort(bootstraps{iModel}.jointlikelihood,'descend');
p = bootstraps{iModel};
parameterEstimates.u0 = p.u0(:,mostLikelyIndices(1));
parameterEstimates.v0 = p.v0(:,mostLikelyIndices(1));
parameterEstimates.u1 = p.u1(:,mostLikelyIndices(1));
parameterEstimates.v1 = p.v1(:,mostLikelyIndices(1));
parameterEstimates.sigma_n = p.sigma_n(:,mostLikelyIndices(1));
parameterEstimates.sigma_s = p.sigma_s(:,mostLikelyIndices(1));
parameterEstimates.zeta = p.zeta(:,mostLikelyIndices(1));
parameterEstimates.delta = p.delta(:,mostLikelyIndices(1));

[u_meso,v_meso,u_bg,v_bg,u_sm,v_sm,dmxdt,dmydt] = DecomposeTrajectories(x, y, t, parameterEstimates);

cv_bg = u_bg + sqrt(-1)*v_bg;
cv_strain = u_meso + sqrt(-1)*v_meso;
cv_sm = u_sm + sqrt(-1)*v_sm;

x_meso = x(1,:) + cumtrapz(t,u_meso);
y_meso = y(1,:) + cumtrapz(t,v_meso);
x_bg = cumtrapz(t,u_bg);
y_bg = cumtrapz(t,v_bg);
x_sm = cumtrapz(t,u_sm);
y_sm = cumtrapz(t,v_sm);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Flyover reference frame
%
% The idea is to visualize the drifters as if we were in a helicoptor
% flying directly above the drifters. So, we will center the window at the
% mesoscale center-of-mass. The window will thus not follow the inertial
% oscillations. 
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x_camera = mean(x_meso,2);
y_camera = mean(y_meso,2);

dx = x - x_camera;
dy = y - y_camera;

buffer = 1.1;
maxDx = buffer*max(dx(:));
minDx = buffer*min(dx(:));
maxDy = buffer*max(dy(:));
minDy = buffer*min(dy(:));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesoscale, fixed coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = figure('Units', 'points', 'Position', [50 50 1920/2 1080/2]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

for iTime=1:1:length(t)
% for iTime=200
    clf 
    
tailLength = 12; % 12 is 6 hours
startIndex = iTime-tailLength;
if (startIndex < 1) startIndex=1; end
range=startIndex:iTime;

s = 1000;
% sp1 = subplot(2,2,[1 3]);
% plot(x/s,y/s,'Color',0.3*[1 1 1]), hold on
sp1 = subplot(1,2,1);
plot(x/s,y/s), axis equal, hold on
rectangle('Position',[x_camera(iTime)+minDx y_camera(iTime)+minDy maxDx-minDx maxDy-minDy]/s,'LineWidth',2)
scatter( (x(iTime,:))/s, (y(iTime,:))/s, 3^2, 'k', 'fill')
sp1.XGrid = 'on';
sp1.YGrid = 'on';
sp1.XMinorGrid = 'on';
sp1.YMinorGrid = 'on';
xlim([-2.5 20])
ylim([-2.5 37.5])

sp2=subplot(1,2,2);
if iTime > 1
    plot( (x(range,:)-0*x_camera(iTime))/s, (y(range,:)-0*y_camera(iTime))/s, 'LineWidth', 2)
    hold on
end
scatter( (x(iTime,:)-0*x_camera(iTime))/s, (y(iTime,:)-0*y_camera(iTime))/s, 8^2, 'k', 'fill')
axis equal
xlim((x_camera(iTime)+[minDx maxDx])/s)
ylim((y_camera(iTime)+[minDy maxDy])/s)
sp2.XGrid = 'on';
sp2.YGrid = 'on';
sp2.XMinorGrid = 'on';
sp2.YMinorGrid = 'on';
sp2.GridAlpha = 0.5;
sp2.MinorGridAlpha = 0.5;
sp2.XTick = -5:5:25;
sp2.YTick = 0:5:35;
sp2.XTickLabel = [];
sp2.YTickLabel = [];
title( sprintf('LatMix 2011 Site 1 Day %d at %2d:%02d', floor(t(iTime)/86400), floor(mod(t(iTime)/3600,24)), floor(mod(t(iTime)/60,60)) ), 'fontsize', 24, 'FontName', 'Helvetica' );

b = 0.05;
p1 = sp1.Position;
sp1.Position = [b/2 b p1(3) 1-2*b];
p1 = sp1.Position;
p2 = sp2.Position;
sp2.Position = [p1(1)+p1(3) p1(2) 1-p1(3)-b 1-2*b];

a = annotation('textarrow',[p1(1)+p1(3) b],[0.075 b],'string',{' Separating Mesoscale and Submesoscale Flows from Clustered Drifter Trajectories.','Oscroft, Sykulski, & Early (2021)'}, 'HeadStyle','none','LineStyle', 'none', 'TextRotation',0);
a.FontSize = 16;
a.FontName = 'Times';

if shouldSaveFigures == 1
    output = sprintf('%s/t_%03d', framesFolder,iTime-1);
%     print('-depsc2', output)
print('-dpng','-r192', output)
end

end
