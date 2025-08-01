% Load our giant list of possible paths.
SiteNumber=1;
K = 4;
dof = 4;
totalPermutations = 1000;
shouldSaveFigures = 0;

load(sprintf('smoothedGriddedRho%dDrifters.mat',SiteNumber));
load(sprintf('BootstrapData/Rho%dDrifterSplineFits%d_dof%d.mat',SiteNumber,totalPermutations,dof));

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

cmp = colormap(parula(7));

errorColor = 0.0*[1 1 1];
errorAlpha = 0.2;
% error68Color = 0.0*[1 1 1];
% error68Alpha = 0.2;
timeMeanColor = 0.2*[1 1 1];
stdEdgeColor = 'none';
meanColor = [0    0.4470    0.7410];
fixedColor = [0.8500    0.3250    0.0980];
meanLineWidth = 2.0;

iModel = 3;

% Values from COM estimates, fit to the entire dataset
sigmaFixed = 0.0591*f0;
thetaFixed = -27.8;

[~,mostLikelyIndices] = sort(bootstraps{iModel}.jointlikelihood,'descend');
sigma = sqrt(bootstraps{iModel}.sigma_n.^2 + bootstraps{iModel}.sigma_s.^2);
theta = atan2(bootstraps{iModel}.sigma_s,bootstraps{iModel}.sigma_n)/2;

indices90 = mostLikelyIndices(1:round(.9*totalPermutations));
indices68 = mostLikelyIndices(1:round(.68*totalPermutations));
indexBest = mostLikelyIndices(1);

error90 = @(value) [max(value(:,indices90),[],2); flip(min(value(:,indices90),[],2))];
error68 = @(value) [max(value(:,indices68),[],2); flip(min(value(:,indices68),[],2))];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Parameter estimate figure
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure('Units', 'points', 'Position', [50 50 figure_width_medium 450*scaleFactor])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

sp4 = subplot(4,1,4);
plot(t/86400,B,'LineWidth',2);
xlim([min(t) max(t)]/86400)
% ylabel('b-splines', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
sp4.YTick = [];
sp4.FontSize = figure_axis_tick_size;
xlabel('time (days)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
c = gray(8);
sp4.ColorOrder = c(4:5,:);

sp0 = subplot(4,1,3);
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
ylim([-0.25 0.25])
lgd = legend([u0p,v0p],'u_0','v_0','Location','Southwest');
lgd.NumColumns = 2;
sp0.FontSize = figure_axis_tick_size;
p = sp0.YLabel.Position;
sp0.YLabel.Position = [p(1)+0.1 p(2) p(3)];

sp1 = subplot(4,1,1);
f = fill([t; flip(t)]/86400,error90(sigma)/f0,errorColor,'EdgeColor',stdEdgeColor); hold on
f.FaceAlpha = errorAlpha;
f = fill([t; flip(t)]/86400,error68(sigma)/f0,errorColor,'EdgeColor',stdEdgeColor);
f.FaceAlpha = errorAlpha;
plot([min(t) max(t)]/86400,sigmaFixed*[1 1]/f0,'k--', 'LineWidth', 1.0*scaleFactor)
sp1.ColorOrderIndex = 1;
plot(t/86400,sigma(:,indexBest)/f0, 'LineWidth', meanLineWidth*scaleFactor)
ylabel('\sigma (f_0)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlim([min(t) max(t)]/86400)
sp1.FontSize = figure_axis_tick_size;
p = sp1.YLabel.Position;
sp1.YLabel.Position = [p(1)+0.1 p(2) p(3)];

sp2 = subplot(4,1,2);
f=fill([t; flip(t)]/86400,error90(theta)*180/pi,errorColor,'EdgeColor',stdEdgeColor); hold on
f.FaceAlpha = errorAlpha;
f=fill([t; flip(t)]/86400,error68(theta)*180/pi,errorColor,'EdgeColor',stdEdgeColor);
f.FaceAlpha = errorAlpha;
plot([min(t) max(t)]/86400,thetaFixed*[1 1],'k--', 'LineWidth', 1.0*scaleFactor)
sp2.ColorOrderIndex = 1;
plot(t/86400,theta(:,indexBest)*180/pi, 'LineWidth', meanLineWidth*scaleFactor)
ylabel('\theta (°)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlim([min(t) max(t)]/86400)
ylim([-100 100])
yticks([-90 -45 0 45 90])
sp2.FontSize = figure_axis_tick_size;


packfig(4,1)


sp0.XLabel.FontSize = figure_axis_label_size;
sp0.YLabel.FontSize = figure_axis_label_size;
sp1.XLabel.FontSize = figure_axis_label_size;
sp1.YLabel.FontSize = figure_axis_label_size;
sp2.XLabel.FontSize = figure_axis_label_size;
sp2.YLabel.FontSize = figure_axis_label_size;
sp4.XLabel.FontSize = figure_axis_label_size;
sp4.YLabel.FontSize = figure_axis_label_size;

p = sp4.Position;
shrinkage = 0.4;
sp4.Position = [p(1) p(2)+0.4*p(4) p(3) 0.6*p(4)];

tightfig

if shouldSaveFigures == 1
    print('Site1Parameters.eps','-depsc2');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Decomposition figures
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = bootstraps{iModel};
u0 = p.u0(:,indexBest);
v0 = p.v0(:,indexBest);
sigma_n = p.sigma_n(:,indexBest);
sigma_s = p.sigma_s(:,indexBest);
zeta = p.zeta(:,indexBest);
delta = p.delta(:,indexBest);

u_meso = u0 + 0.5*(sigma_n+delta).*x + 0.5*(sigma_s-zeta).*y;
v_meso = v0 + 0.5*(sigma_s+zeta).*x + 0.5*(delta-sigma_n).*y;
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Quick aside
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

measuredTotalVariance = sum(dqdt.^2 + drdt.^2,2);
modelMesoscaleVariance = sum(dqdt_meso.^2 + drdt_meso.^2,2);
modelSubmesoscaleVariance = sum(u_sm.^2 + v_sm.^2,2);
modelTotalVariance = sum(dqdt_meso.^2 + drdt_meso.^2+u_sm.^2 + v_sm.^2,2);
EnergyCrossTerms = (measuredTotalVariance - modelTotalVariance)./measuredTotalVariance;
figure
plot(t/86400,EnergyCrossTerms) %, hold on, plot(t/86400,actualTotalVariance);

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
text(14,33,'meso', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlabel('km', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('km', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
sp1.FontSize = figure_axis_tick_size;
sp1.XLabel.FontSize = figure_axis_label_size;
sp1.YLabel.FontSize = figure_axis_label_size;

if shouldSaveFigures == 1
    print('Site1DecompFigA.eps','-depsc2');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Background and submesoscale
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f2 = figure('Units', 'points', 'Position', [50 50 figure_width_1col figureHeight*scaleFactor]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

s = 1000;

sp2 = subplot(2,1,1);
plot(x_bg/s,y_bg/s), axis equal
text(.75,0.65,'bg', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
sp2.YAxisLocation = 'right';
sp2.XTickLabel = [];
ylabel('km', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
sp2.FontSize = figure_axis_tick_size;

sp3 = subplot(2,1,2);
plot(x_sm/s,y_sm/s), axis equal
text(.75,0.65,'sm', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
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
    print('Site1DecompFigB.eps','-depsc2');
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Mesoscale, COM coordinates
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

f3 = figure('Units', 'points', 'Position', [50 50 figure_width_2col figureHeight*scaleFactor]);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

s = 1000;
% sp1 = subplot(2,2,[1 3]);
% plot(x/s,y/s,'Color',0.3*[1 1 1]), hold on
plot(q_meso/s,r_meso/s,'LineWidth',1.5), axis equal
text(0.75,1.75,'meso (centre-of-mass)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
xlabel('km', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('km', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
sp1.FontSize = figure_axis_tick_size;
sp1.XLabel.FontSize = figure_axis_label_size;
sp1.YLabel.FontSize = figure_axis_label_size;

if shouldSaveFigures == 1
    print('Site1DecompFigC.eps','-depsc2');
end

