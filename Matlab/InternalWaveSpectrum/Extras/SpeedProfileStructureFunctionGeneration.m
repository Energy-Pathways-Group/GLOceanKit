profilename = 'latmix-site1';
% profilename = 'exponential';
boundaryCondition = UpperBoundary.rigidLid;
latitude = 31;

[rho, N2, zIn] = InternalModes.StratificationProfileWithName(profilename);

GM = GarrettMunkSpectrum(rho,zIn,latitude);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute and plot the variances as a function of depth.
%
z = linspace(min(zIn),max(zIn),5000)';

approximation = 'gm';
approximation = 'exact';
Euv = GM.HorizontalVelocityVariance(z,approximation);
Eeta = GM.IsopycnalVariance(z,approximation);
Ew = GM.VerticalVelocityVariance(z,approximation);
E = 0.5*(Euv + Ew + N2(z).*Eeta); % total energy

fprintf('Total depth integrated energy: %f\n',trapz(z,E));

figure
subplot(1,4,1)
plot(1e4*Euv,z)
xlabel('velocity variance (cm^2/s^2)')
ylabel('depth (m)')
subplot(1,4,2)
plot(Eeta,z)
xlabel('isopycnal variance (m^2)')
set(gca,'YTickLabel',[]);
subplot(1,4,3)
plot(1e4*Ew,z)
xlabel('vertical velocity variance (cm^2/s^2)')
set(gca,'YTickLabel',[]);
subplot(1,4,4)
plot(1e4*E,z)
xlabel('total energy (cm^2/s^2)')
set(gca,'YTickLabel',[]);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Vertical structure functions
%

scaleFactor = 1;
LoadFigureDefaults;


N0 = sqrt(max(N2(z)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% pretty axis labels
%



shouldUseLogForOmega = 1;

xAxisMax = N0;
f0 = GM.f0;

if shouldUseLogForOmega == 1
    ticks_x = exp(linspace(log(f0),log(xAxisMax),5));
else
    ticks_x = linspace(f0,xAxisMax,5);
end

labels_x = cell(length(ticks_x),1);
for i=1:length(ticks_x)
    if round(ticks_x(i)/f0) == 0.0
        labels_x{i} = '0';
    elseif round(ticks_x(i)/f0) == 1
        labels_x{i} = 'f_0';
    elseif round(ticks_x(i)/f0) == -1
        labels_x{i} = '-f_0';
    elseif ticks_x(i) == N0
        labels_x{i} = 'N_0';
    elseif ticks_x(i) == -N0
        labels_x{i} = '-N_0';
    else
        labels_x{i} = sprintf('%df_0',round(ticks_x(i)/f0));
    end
end

%%%% Small fudge factor, move the tick, keep the label name.
ticks_x(end) = max(GM.omega);

% Now for depth
shouldUseLogForDepth = 1;
if shouldUseLogForDepth == 1
    axis_depth = -log(abs(GM.zInternal)+10);
    ticks_y = linspace(min(axis_depth),max(axis_depth),6);
    ticks_y_depths = round(-(exp(-ticks_y)-10));
else
    axis_depth = GM.zInternal;
    ticks_y = linspace(min(axis_depth),max(axis_depth),6);
    ticks_y_depths = ticks_y;
end
labels_y = cell(length(ticks_y),1);
for i=1:length(ticks_y)
    labels_y{i} = sprintf('%d',ticks_y_depths(i));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Turning point function
%

D = max(zIn)-min(zIn);
b = GM.L_gm;
N_gm = GM.invT_gm;

N = sqrt(N2(GM.zInternal));
Phi = GM.Phi_omega;
Gamma = GM.Gamma_omega;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Vertical structure function figure
%
contour_color = 0.5*[1 1 1];
contour_width = 1.0;

% change this function between log10(z) and (z) to change the scale
f = @(z) log10(z);

% wkb normalized version
% PhiPlot = b*N_gm*Phi./N;
% GammaPlot = b*N_gm*Gamma.*N;

% unnormalized version
PhiPlot = sum(Phi,3,'omitnan');
GammaPlot = sum(Gamma,3,'omitnan');

% Now apply the function
PhiVariance = f(PhiPlot);
GammaVariance = f(GammaPlot);

% We the first column (omega=f_0) to normalize.
PhiVariance_range = f( [1e-4 1]*max(PhiPlot(:,1)) );
GammaVariance_range = f( [1e-4 1]*max(GammaPlot(:,1)) );

FigureSize = [50 50 figure_width_1col+8 225*scaleFactor];
fig1 = figure('Units', 'points', 'Position', FigureSize);
% set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

p1 = subplot(1,2,1);
pcolor( GM.omega, axis_depth, PhiVariance ), xlog, shading flat, hold on
plot( sqrt(N2(GM.zInternal)), axis_depth, 'LineWidth', 2.0*scaleFactor, 'Color', 1.0*[1 1 1])
caxis(PhiVariance_range)
set( gca, 'FontSize', figure_axis_tick_size);
set(gca, 'YTick', 1000*(-5:1:0));
ylabel('depth (m)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
title('$\Phi(\omega,z)$','Interpreter','LaTex', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
xticks(ticks_x)
xticklabels(labels_x)
yticks(ticks_y)
yticklabels(labels_y)

[C,~] = contour(GM.omega, axis_depth, PhiVariance,[1 1]*min(PhiVariance_range));
startIndex = 1;
while (startIndex < size(C,2))
    endIndex = startIndex+C(2,startIndex);
    plot( C(1,(startIndex+1):endIndex), C(2,(startIndex+1):endIndex),'LineStyle', '--', 'LineWidth', contour_width*scaleFactor, 'Color', contour_color )
    startIndex = endIndex + 1;
end

p2 = subplot(1,2,2);
pcolor( GM.omega, axis_depth, GammaVariance ),xlog , shading flat, hold on
plot( sqrt(N2(GM.zInternal)), axis_depth, 'LineWidth', 2.0*scaleFactor, 'Color', 1.0*[1 1 1])
caxis(GammaVariance_range)
set( gca, 'FontSize', figure_axis_tick_size);
set(gca, 'YTick', []);
title('$\Gamma(\omega,z)$','Interpreter','LaTex', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
xticks(ticks_x(2:end))
xticklabels(labels_x(2:end))

[C,~] = contour(GM.omega, axis_depth, GammaVariance,[1 1]*min(GammaVariance_range));
startIndex = 1;
while (startIndex < size(C,2))
    endIndex = startIndex+C(2,startIndex);
    plot( C(1,(startIndex+1):endIndex), C(2,(startIndex+1):endIndex),'LineStyle', '--', 'LineWidth', contour_width*scaleFactor, 'Color', contour_color )
    startIndex = endIndex + 1;
end


colormap( cmocean('dense') );

cb = colorbar('eastoutside');
cb.Ticks = f([1e-4 1e-3 1e-2 1e-1 1]*max(GammaPlot(:,1)));
cb.TickLabels = {'10^{-4}', '10^{-3}', '10^{-2}', '10^{-1}', '10^0'};

pause(0.1) % Bizarrely, if you don't pause, the position doesn't work correctly. This is insane.
p1widthIncrease = 0.045;
p1.Position = [p1.Position(1) p1.Position(2) p1.Position(3)+p1widthIncrease p1.Position(4)];
space = 0.02;
p2.Position = [p1.Position(1)+p1.Position(3)+space p1.Position(2) p1.Position(3) p1.Position(4)];
cb.Position = [p2.Position(1)+p2.Position(3)+space cb.Position(2) cb.Position(3) cb.Position(4)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Vertical structure function figure for K
%
Phi = GM.Phi_k;
Gamma = GM.Gamma_k;

% change this function between log10(z) and (z) to change the scale
f = @(z) log10(z);

% unnormalized version
PhiPlot = sum(Phi,3,'omitnan');
GammaPlot = sum(Gamma,3,'omitnan');

% Now apply the function
PhiVariance = f(PhiPlot);
GammaVariance = f(GammaPlot);

PhiVariance_range = f( [1e-4 1]*max(PhiPlot(:,1)) );
GammaVariance_range = f( [1e-4 1]*max(GammaPlot(:,1)) );

FigureSize = [50 50 figure_width_1col+8 225*scaleFactor];
fig1 = figure('Units', 'points', 'Position', FigureSize);
% set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

p1 = subplot(1,2,1);
pcolor( GM.k, axis_depth, PhiVariance ), xlog, shading flat, hold on
caxis(PhiVariance_range)
set( gca, 'FontSize', figure_axis_tick_size);
set(gca, 'YTick', 1000*(-5:1:0));
ylabel('depth (m)', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
title('$\Phi(k,z)$','Interpreter','LaTex', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
% xticks(ticks_x)
%xticklabels(labels_x)
yticks(ticks_y)
yticklabels(labels_y)

p2 = subplot(1,2,2);
pcolor( GM.k, axis_depth, GammaVariance ),xlog , shading flat, hold on
caxis(GammaVariance_range)
set( gca, 'FontSize', figure_axis_tick_size);
set(gca, 'YTick', []);
title('$\Gamma(k,z)$','Interpreter','LaTex', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
% xticks(ticks_x(2:end))
% xticklabels(labels_x(2:end))


return

if boundaryCondition == UpperBoundary.freeSurface
    filename = sprintf('%s-free-surface.mat',profilename);
else
    filename = sprintf('%s.mat',profilename);
end

save(filename,'F_omega','G_omega','omega', 'h_omega', 'k_omega', 'F_k','G_k','h_k', 'k', 'omega_k', 'latitude','zIn','N_max','zInternal','N2internal', 'rho', 'N2');