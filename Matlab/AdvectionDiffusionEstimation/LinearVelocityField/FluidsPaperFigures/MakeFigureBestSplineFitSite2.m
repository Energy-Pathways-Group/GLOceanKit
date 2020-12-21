% Load our giant list of possible paths.
SiteNumber=2;
dof = 1;
totalPermutations = 1000;
shouldSaveFigures = 0;

load(sprintf('smoothedGriddedRho%dDrifters.mat',SiteNumber));
load(sprintf('BootstrapData/Rho%dDrifterSplineFits%d_dof%d-nou1v1.mat',SiteNumber,totalPermutations,dof));
% load(sprintf('BootstrapData/Rho%dDrifterSplineFits%d_dof%d.mat',SiteNumber,totalPermutations,dof));

scaleFactor = 1;
LoadFigureDefaults

f0 = 2 * 7.2921E-5 * sin( lat0*pi/180. );
x = x(:,1:(end-1));
y = y(:,1:(end-1));
nDrifters = size(x,2);
nT = size(x,1);
D = FiniteDifferenceMatrix(1,t,1,1,2);
mx = mean(x,2);
my = mean(y,2);
dmxdt = D*mx;
dmydt = D*my;
dxdt = D*x;
dydt = D*y;


errorColor = 0.0*[1 1 1];
errorAlpha = 0.2;
timeMeanColor = 0.2*[1 1 1];
stdEdgeColor = 'none';
meanColor = [0 0 1];
fixedColor = [1 0 0];
meanLineWidth = 2.0;

iModel = 3;

% Values from COM estimates, fit to the entire dataset
sigmaFixed = 0.0642*f0;
thetaFixed = -78;
zetaFixed = 0.065*f0;

[~,mostLikelyIndices] = sort(bootstraps{iModel}.jointlikelihood,'descend');
sigma = sqrt(bootstraps{iModel}.sigma_n.^2 + bootstraps{iModel}.sigma_s.^2);
theta = atan2(bootstraps{iModel}.sigma_s,bootstraps{iModel}.sigma_n)/2;
zeta = bootstraps{iModel}.zeta;
s2 = sigma.^2-zeta.^2;

indices90 = mostLikelyIndices(1:round(.9*totalPermutations));
indices68 = mostLikelyIndices(1:round(.68*totalPermutations));
indexBest = mostLikelyIndices(1);

error90 = @(value) [max(value(:,indices90),[],2); flip(min(value(:,indices90),[],2))];
error68 = @(value) [max(value(:,indices68),[],2); flip(min(value(:,indices68),[],2))];

% for i=1:totalPermutations
%    if ( 
% end

for i=1:totalPermutations
   flipID = find(abs(diff(theta(:,i)*180/pi)) > 90, 1, 'first');
   theta(1:flipID,i) = theta(1:flipID,i)-pi;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parameter estimate figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Units', 'points', 'Position', [50 50 figure_width_medium 550*scaleFactor])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

sp4 = subplot(5,1,5);
plot(t/86400,B,'LineWidth',2);
xlim([min(t) max(t)]/86400)
% ylabel('b-splines', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
sp4.YTick = [];
sp4.FontSize = figure_axis_tick_size;
xlabel('time (days)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
c = gray(8);
sp4.ColorOrder = c(4:5,:);

% p = bootstraps{iModel};
% u_com_meso = p.u0 + 0.5*(p.sigma_n+p.delta).*mx + 0.5*(p.sigma_s-p.zeta).*my;
% v_com_meso = p.v0 + 0.5*(p.sigma_s+p.zeta).*mx + 0.5*(p.delta-p.sigma_n).*my;

sp0 = subplot(5,1,4);
f = fill([t; flip(t)]/86400,error90(bootstraps{iModel}.u0),errorColor,'EdgeColor',stdEdgeColor); hold on
f.FaceAlpha = errorAlpha;
f = fill([t; flip(t)]/86400,error68(bootstraps{iModel}.u0),errorColor,'EdgeColor',stdEdgeColor); hold on
f.FaceAlpha = errorAlpha;
f = fill([t; flip(t)]/86400,error90(bootstraps{iModel}.v0),errorColor,'EdgeColor',stdEdgeColor);
f.FaceAlpha = errorAlpha;
f = fill([t; flip(t)]/86400,error68(bootstraps{iModel}.v0),errorColor,'EdgeColor',stdEdgeColor);
f.FaceAlpha = errorAlpha;
sp0.ColorOrderIndex = 3;
u0p = plot(t/86400,bootstraps{iModel}.u0(:,indexBest), 'LineWidth', meanLineWidth*scaleFactor);
v0p = plot(t/86400,bootstraps{iModel}.v0(:,indexBest), 'LineWidth', meanLineWidth*scaleFactor);
ylabel('u_0/v_0 (m/s)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlim([min(t) max(t)]/86400)
ylim([-5 5])
lgd = legend([u0p,v0p],'u_0','v_0','Location','southwest');
lgd.NumColumns = 2;
sp0.FontSize = figure_axis_tick_size;
p = sp0.YLabel.Position;
sp0.YLabel.Position = [p(1)+0.1 p(2) p(3)];

s = f0;
sp1 = subplot(5,1,2);
f = fill([t; flip(t)]/86400,error90(sigma)/f0,errorColor,'EdgeColor',stdEdgeColor); hold on
f.FaceAlpha = errorAlpha;
f = fill([t; flip(t)]/86400,error68(sigma)/f0,errorColor,'EdgeColor',stdEdgeColor);
f.FaceAlpha = errorAlpha;
plot([min(t) max(t)]/86400,sigmaFixed*[1 1]/f0,'k--', 'LineWidth', meanLineWidth*scaleFactor)
sp1.ColorOrderIndex = 1;
plot(t/86400,sigma(:,indexBest)/s, 'LineWidth', meanLineWidth*scaleFactor)
ylabel('\sigma (f_0)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlim([min(t) max(t)]/86400)
sp1.FontSize = figure_axis_tick_size;
ylim([0 0.8])
% p = sp1.YLabel.Position;
% sp1.YLabel.Position = [p(1)+0.2 p(2) p(3)];

sp2 = subplot(5,1,3);
f=fill([t; flip(t)]/86400,error90(theta)*180/pi,errorColor,'EdgeColor',stdEdgeColor); hold on
f.FaceAlpha = errorAlpha;
f=fill([t; flip(t)]/86400,error68(theta)*180/pi,errorColor,'EdgeColor',stdEdgeColor);
f.FaceAlpha = errorAlpha;
plot([min(t) max(t)]/86400,thetaFixed*[1 1],'k--', 'LineWidth', meanLineWidth*scaleFactor)
sp2.ColorOrderIndex = 1;
plot(t/86400,theta(:,indexBest)*180/pi, 'LineWidth', meanLineWidth*scaleFactor)
ylabel('\theta (°)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlim([min(t) max(t)]/86400)
ylim([-140 40])
yticks([-90 -45 0 45 90])
sp2.FontSize = figure_axis_tick_size;
xlabel('time (days)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)

sp3 = subplot(5,1,1);
f = fill([t; flip(t)]/86400,error90(zeta)/f0,errorColor,'EdgeColor',stdEdgeColor); hold on
f.FaceAlpha = errorAlpha;
f = fill([t; flip(t)]/86400,error68(zeta)/f0,errorColor,'EdgeColor',stdEdgeColor);
f.FaceAlpha = errorAlpha;
plot([min(t) max(t)]/86400,zetaFixed*[1 1],'k--', 'LineWidth', meanLineWidth*scaleFactor)
sp3.ColorOrderIndex = 7;
plot(t/86400,zeta(:,indexBest)/f0, 'LineWidth', meanLineWidth*scaleFactor)
ylabel('\zeta (f_0)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlim([min(t) max(t)]/86400)
ylim([-1 1])
sp3.FontSize = figure_axis_tick_size;
xlabel('time (days)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)

packfig(5,1)

sp0.XLabel.FontSize = figure_axis_label_size;
sp0.YLabel.FontSize = figure_axis_label_size;
sp1.XLabel.FontSize = figure_axis_label_size;
sp1.YLabel.FontSize = figure_axis_label_size;
sp2.XLabel.FontSize = figure_axis_label_size;
sp2.YLabel.FontSize = figure_axis_label_size;
sp3.XLabel.FontSize = figure_axis_label_size;
sp3.YLabel.FontSize = figure_axis_label_size;
sp4.XLabel.FontSize = figure_axis_label_size;
sp4.YLabel.FontSize = figure_axis_label_size;

p = sp4.Position;
shrinkage = 0.4;
sp4.Position = [p(1) p(2)+0.4*p(4) p(3) 0.6*p(4)];

tightfig

if shouldSaveFigures == 1
    print('Site2Parameters.eps','-depsc2');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parameter estimate figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = bootstraps{iModel};
u0 = p.u0(:,indexBest);
v0 = p.v0(:,indexBest);
u1 = 0; %p.u1(:,indexBest);
v1 = 0; %p.v1(:,indexBest);
sigma_n = p.sigma_n(:,indexBest);
sigma_s = p.sigma_s(:,indexBest);
zeta = p.zeta(:,indexBest);
delta = p.delta(:,indexBest);

u_meso = u0 + u1.*t + 0.5*(sigma_n+delta).*x + 0.5*(sigma_s-zeta).*y;
v_meso = v0 + v1.*t + 0.5*(sigma_s+zeta).*x + 0.5*(delta-sigma_n).*y;
x_meso = x(1,:) + cumtrapz(t,u_meso);
y_meso = y(1,:) + cumtrapz(t,v_meso);

u_bg = dmxdt - mean(u_meso,2);
v_bg = dmydt - mean(v_meso,2);
x_bg = cumtrapz(t,u_bg);
y_bg = cumtrapz(t,v_bg);

u_sm = dxdt - u_meso - u_bg;
v_sm = dydt - v_meso - v_bg;
x_sm = cumtrapz(t,u_sm);
y_sm = cumtrapz(t,v_sm);

q = (x-mx);
r = (y-my);
dqdt = D*q;
drdt = D*r;
dqdt_meso = 0.5*(sigma_n+delta).*q + 0.5*(sigma_s-zeta).*r;
drdt_meso = 0.5*(sigma_s+zeta).*q + 0.5*(delta-sigma_n).*r;
q_meso = q(1,:) + cumtrapz(t,dqdt_meso);
r_meso = r(1,:) + cumtrapz(t,drdt_meso);

figureHeight = 367;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesoscale, fixed coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = figure('Units', 'points', 'Position', [50 50 figure_width_1col figureHeight*scaleFactor]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

s = 1000;
% sp1 = subplot(2,2,[1 3]);
% plot(x/s,y/s,'Color',0.3*[1 1 1]), hold on
plot(x_meso/s,y_meso/s), axis equal
text(75,200,'meso', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlabel('km', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('km', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
sp1.FontSize = figure_axis_tick_size;
sp1.XLabel.FontSize = figure_axis_label_size;
sp1.YLabel.FontSize = figure_axis_label_size;

if shouldSaveFigures == 1
    print('Site2DecompFigA.eps','-depsc2');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesoscale, COM coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f1 = figure('Units', 'points', 'Position', [50 50 figure_width_1col figureHeight*scaleFactor]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

s = 1000;
% sp1 = subplot(2,2,[1 3]);
% plot(x/s,y/s,'Color',0.3*[1 1 1]), hold on
plot(q_meso/s,r_meso/s,'LineWidth',1.5), axis equal
text(-0.5,4,'meso (centre-of-mass)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlabel('km', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('km', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
sp1.FontSize = figure_axis_tick_size;
sp1.XLabel.FontSize = figure_axis_label_size;
sp1.YLabel.FontSize = figure_axis_label_size;

if shouldSaveFigures == 1
    print('Site2DecompFigC.eps','-depsc2');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Background and submesoscale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f2 = figure('Units', 'points', 'Position', [50 50 figure_width_1col figureHeight*scaleFactor]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

s = 1000;

xtext = 2;
ytext = 3;

sp2 = subplot(2,1,1);
plot(x_bg/s,y_bg/s), axis equal
text(xtext,ytext,'bg', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
sp2.YAxisLocation = 'right';
sp2.XTickLabel = [];
ylabel('km', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
sp2.FontSize = figure_axis_tick_size;

sp3 = subplot(2,1,2);
plot(x_sm/s,y_sm/s), axis equal
text(xtext,ytext,'sm', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
sp3.YAxisLocation = 'right';
xlabel('km', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('km', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
sp3.FontSize = figure_axis_tick_size;

sp2sp3xlim = [min([x_bg(:);x_sm(:)]) max([x_bg(:);x_sm(:)])]/s;
range = max(sp2sp3xlim)-min(sp2sp3xlim);
sp2sp3xlim = [sp2sp3xlim(1)-0.05*range sp2sp3xlim(2)+0.05*range];

sp2.XLim = sp2sp3xlim;
sp3.XLim = sp2sp3xlim;

sp2sp3ylim = [min([y_bg(:);y_sm(:)]) max([y_bg(:);y_sm(:)])]/s;
range = max(sp2sp3ylim)-min(sp2sp3ylim);
sp2sp3ylim = [sp2sp3ylim(1)-0.05*range sp2sp3ylim(2)+0.05*range];

sp2.YLim = sp2sp3ylim;
sp3.YLim = sp2sp3ylim;

packfig(2,1)

sp2.XLabel.FontSize = figure_axis_label_size;
sp2.YLabel.FontSize = figure_axis_label_size;
sp3.XLabel.FontSize = figure_axis_label_size;
sp3.YLabel.FontSize = figure_axis_label_size;

if shouldSaveFigures == 1
    print('Site2DecompFigB.eps','-depsc2');
end

