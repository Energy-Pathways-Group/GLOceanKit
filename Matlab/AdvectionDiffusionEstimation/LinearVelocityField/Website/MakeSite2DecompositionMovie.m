% Load our giant list of possible paths.
SiteNumber=2;
dof = 6;
iModel = 3;
totalPermutations = 1000;
shouldSaveFigures = 1;
framesFolder = 'Site2DecompositionFramesModel3Dof6';

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

buffer = 1.0;
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
% for iTime=150
    clf 
    
tailLength = 12; % 12 is 6 hours
startIndex = iTime-tailLength;
if (startIndex < 1) startIndex=1; end
range=startIndex:iTime;

s = 1000;

sp1=subplot(6,3,[2 5 8 11 14 17]);
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
sp1.XTick = -25:5:125;
sp1.YTick = 0:5:230;
sp1.XTickLabel = [];
sp1.YTickLabel = [];
title( sprintf('LatMix 2011 Site 2 Day %d at %2d:%02d', floor(t(iTime)/86400), floor(mod(t(iTime)/3600,24)), floor(mod(t(iTime)/60,60)) ), 'fontsize', 24, 'FontName', 'Helvetica' );


sp2 = subplot(6,3,[3 6 9]);
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
sp2.Box = 'on';
text(-2.75,2.5,'background', 'FontSize', 16, 'FontName', 'Helvetica')

sp3 = subplot(6,3,[12 15 18]);
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
sp3.Box = 'on';
text(-2.75,2.5,'submesoscale', 'FontSize', 16, 'FontName', 'Helvetica')

sp2.XLim = 3.2*[-1 1];
sp3.XLim = 3.2*[-1 1];

sp2.YLim = 3.2*[-1 1];
sp3.YLim = 3.2*[-1 1];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Parameter estimates and plots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
f0 = 2 * 7.2921E-5 * sin( lat0*pi/180. );

errorColor = 0.0*[1 1 1];
errorAlpha = 0.2;
timeMeanColor = 0.2*[1 1 1];
stdEdgeColor = 'none';
meanColor = [0 0 1];
fixedColor = [1 0 0];
meanLineWidth = 2.0;
scaleFactor = 1.5;
LoadFigureDefaults

sigma = sqrt(bootstraps{iModel}.sigma_n.^2 + bootstraps{iModel}.sigma_s.^2);
theta = atan2(bootstraps{iModel}.sigma_s,bootstraps{iModel}.sigma_n)/2;
zeta = bootstraps{iModel}.zeta;
s2 = sigma.^2-zeta.^2;

indices90 = mostLikelyIndices(1:round(.9*totalPermutations));
indices68 = mostLikelyIndices(1:round(.68*totalPermutations));
indexBest = mostLikelyIndices(1);

error90 = @(value) [max(value(:,indices90),[],2); flip(min(value(:,indices90),[],2))];
error68 = @(value) [max(value(:,indices68),[],2); flip(min(value(:,indices68),[],2))];

for i=1:totalPermutations
   flipID = find(abs(diff(theta(:,i)*180/pi)) > 90, 1, 'first');
   theta(1:flipID,i) = theta(1:flipID,i)-pi;
end

s = f0;
sp4 = subplot(6,3,[1 4]);
f = fill([t; flip(t)]/86400,error90(sigma)/f0,errorColor,'EdgeColor',stdEdgeColor); hold on
f.FaceAlpha = errorAlpha;
f = fill([t; flip(t)]/86400,error68(sigma)/f0,errorColor,'EdgeColor',stdEdgeColor);
f.FaceAlpha = errorAlpha;
sp4.ColorOrderIndex = 1;
plot(t/86400,sigma(:,indexBest)/s, 'LineWidth', meanLineWidth*scaleFactor)
ylabel('\sigma (f_0)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlim([min(t) max(t)]/86400)
sp4.FontSize = figure_axis_tick_size;
sp4.XTickLabel = [];
ylimits = [0 0.8];
ylim(ylimits)
plot([t(iTime) t(iTime)]/86400, ylimits,'LineWidth',2,'Color',0*[1 1 1])
% p = sp1.YLabel.Position;
% sp1.YLabel.Position = [p(1)+0.2 p(2) p(3)];

sp5 = subplot(6,3,[7 10]);
f=fill([t; flip(t)]/86400,error90(theta)*180/pi,errorColor,'EdgeColor',stdEdgeColor); hold on
f.FaceAlpha = errorAlpha;
f=fill([t; flip(t)]/86400,error68(theta)*180/pi,errorColor,'EdgeColor',stdEdgeColor);
f.FaceAlpha = errorAlpha;
sp5.ColorOrderIndex = 1;
plot(t/86400,theta(:,indexBest)*180/pi, 'LineWidth', meanLineWidth*scaleFactor)
ylabel('\theta (°)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlim([min(t) max(t)]/86400)
ylimits = [-140 40];
ylim(ylimits)
plot([t(iTime) t(iTime)]/86400, ylimits,'LineWidth',2,'Color',0*[1 1 1])
yticks([-90 -45 0 45 90])
sp5.FontSize = figure_axis_tick_size;
% xlabel('time (days)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
sp5.XTickLabel = [];

sp6 = subplot(6,3,[13 16]);
f = fill([t; flip(t)]/86400,error90(zeta)/f0,errorColor,'EdgeColor',stdEdgeColor); hold on
f.FaceAlpha = errorAlpha;
f = fill([t; flip(t)]/86400,error68(zeta)/f0,errorColor,'EdgeColor',stdEdgeColor);
f.FaceAlpha = errorAlpha;
sp6.ColorOrderIndex = 7;
plot(t/86400,zeta(:,indexBest)/f0, 'LineWidth', meanLineWidth*scaleFactor)
ylabel('\zeta (f_0)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlim([min(t) max(t)]/86400)
ylimits = [-1 1];
ylim(ylimits);
plot([t(iTime) t(iTime)]/86400, ylimits,'LineWidth',2,'Color',0*[1 1 1])
sp6.FontSize = figure_axis_tick_size;
xlabel('time (days)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)



b = 0.07;
dx = 0.06;
dy = -0.005;
w=0.3;
h=0.27;
sp4.Position = [dx b+2*dy+2*0.3 w h];
sp5.Position = [dx b+dy+0.3 w h];
sp6.Position = [dx b w h];





sp1.Position = [0.3 b 0.5 1-2*b];
p1 = sp1.Position;
p2 = sp2.Position;
w=0.42;
h=w;
sp2.Position = [0.64 0.51 w h];
sp3.Position = [0.64 0.07 w h];
% b = 0.05;
% p1 = sp1.Position;
% sp1.Position = [b/2 b p1(3) 1-2*b];
% p1 = sp1.Position;
% p2 = sp2.Position;
% sp2.Position = [p1(1)+p1(3) p1(2) 1-p1(3)-b 1-2*b];
% 
a = annotation('textarrow',0.97*[1 1],0.04*[1 1],'string',{'Separating Mesoscale and Submesoscale Flows from Clustered Drifter Trajectories.','Oscroft, Sykulski, & Early (2021)'}, 'HeadStyle','none','LineStyle', 'none', 'TextRotation',0);
a.FontSize = 16;
a.FontName = 'Times';

b = annotation('textarrow',0.48*[1 1],0.92*[1 1],'string','mesoscale', 'HeadStyle','none','LineStyle', 'none', 'TextRotation',0);
b.FontSize = 16;
b.FontName = 'Helvetica';

if shouldSaveFigures == 1
    output = sprintf('%s/t_%03d', framesFolder,iTime-1);
%     print('-depsc2', output)
print('-dpng','-r192', output)
end

end
