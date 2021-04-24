% Load our giant list of possible paths.
SiteNumber=1;
K = 4;
dof = 4;
iModel = 3;
totalPermutations = 1000;
shouldSaveFigures = 1;
framesFolder = 'Site1DecompositionFrames';

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

% for iTime=1:1:length(t)
for iTime=1
    clf 
    
tailLength = 12; % 12 is 6 hours
startIndex = iTime-tailLength;
if (startIndex < 1) startIndex=1; end
range=startIndex:iTime;

s = 1000;

sp1=subplot(2,3,[1 2 4 5]);
plot( x/s, y/s, 'LineWidth', 1, 'Color', 0.8*[1 1 1])
hold on
if iTime > 1
    plot( x_meso(range,:)/s, y_meso(range,:)/s, 'LineWidth', 2)
end
% scatter( x(iTime,:)/s, y(iTime,:)/s, 8^2, 0.4*[1 1 1], 'fill'), hold on
scatter( x_meso(iTime,:)/s, y_meso(iTime,:)/s, 8^2, 'k', 'fill')

axis equal
xlim((x_camera(iTime)+[minDx maxDx])/s)
ylim((y_camera(iTime)+[minDy maxDy])/s)
sp1.XGrid = 'on';
sp1.YGrid = 'on';
sp1.XMinorGrid = 'on';
sp1.YMinorGrid = 'on';
sp1.GridAlpha = 0.5;
sp1.MinorGridAlpha = 0.5;
sp1.XTick = -5:5:25;
sp1.YTick = 0:5:35;
sp1.XTickLabel = [];
sp1.YTickLabel = [];
title( sprintf('LatMix 2011 Site 1 Day %d at %2d:%02d', floor(t(iTime)/86400), floor(mod(t(iTime)/3600,24)), floor(mod(t(iTime)/60,60)) ), 'fontsize', 24, 'FontName', 'Helvetica' );


sp2 = subplot(2,3,3);
if iTime > 1
plot(x_bg(range)/s,y_bg(range)/s, 'LineWidth', 2, 'color', 0.4*[1 1 1]), hold on
end
scatter( x_bg(iTime)/s,y_bg(iTime)/s, 8^2, 'k', 'fill'), axis equal
sp2.XGrid = 'on';
sp2.YGrid = 'on';
sp2.XMinorGrid = 'on';
sp2.YMinorGrid = 'on';
sp2.GridAlpha = 0.5;
sp2.MinorGridAlpha = 0.5;
sp2.XTick = -5:5:5;
sp2.YTick = -5:5:5;
sp2.XTickLabel = [];
sp2.YTickLabel = [];
sp2.Box = 'on'
text(-1,1.1,'background', 'FontSize', 16, 'FontName', 'Helvetica')

sp3 = subplot(2,3,6);
if iTime > 1
    plot( x_sm(range,:)/s, y_sm(range,:)/s, 'LineWidth', 2), hold on
end
% scatter( x(iTime,:)/s, y(iTime,:)/s, 8^2, 0.4*[1 1 1], 'fill'), hold on
scatter( x_sm(iTime,:)/s, y_sm(iTime,:)/s, 8^2, 'k', 'fill'), axis equal
sp3.XGrid = 'on';
sp3.YGrid = 'on';
sp3.XMinorGrid = 'on';
sp3.YMinorGrid = 'on';
sp3.GridAlpha = 0.5;
sp3.MinorGridAlpha = 0.5;
sp3.XTick = -5:5:5;
sp3.YTick = -5:5:5;
sp3.XTickLabel = [];
sp3.YTickLabel = [];
sp3.Box = 'on'
text(-1,1.1,'submesoscale', 'FontSize', 16, 'FontName', 'Helvetica')

sp2.XLim = [-1.25 1.25];
sp3.XLim = [-1.25 1.25];

sp2.YLim = [-1.25 1.25];
sp3.YLim = [-1.25 1.25];

b = 0.05;
sp1.Position = [b/2 b 0.66 1-2*b];
p1 = sp1.Position;
p2 = sp2.Position;
sp2.Position = [0.62 0.51 0.41 0.41];
sp3.Position = [0.62 0.07 0.41 0.41];
% b = 0.05;
% p1 = sp1.Position;
% sp1.Position = [b/2 b p1(3) 1-2*b];
% p1 = sp1.Position;
% p2 = sp2.Position;
% sp2.Position = [p1(1)+p1(3) p1(2) 1-p1(3)-b 1-2*b];
% 
a = annotation('textarrow',0.814*[1 1],0.045*[1 1],'string',{' Separating Mesoscale and Submesoscale Flows from Clustered Drifter Trajectories. Oscroft, Sykulski, & Early (2021)'}, 'HeadStyle','none','LineStyle', 'none', 'TextRotation',0);
a.FontSize = 16;
a.FontName = 'Times';

b = annotation('textarrow',0.13*[1 1],0.9*[1 1],'string','mesoscale', 'HeadStyle','none','LineStyle', 'none', 'TextRotation',0);
b.FontSize = 16;
b.FontName = 'Helvetica';

if shouldSaveFigures == 1
    output = sprintf('%s/t_%03d', framesFolder,iTime-1);
%     print('-depsc2', output)
print('-dpng','-r192', output)
end

end
