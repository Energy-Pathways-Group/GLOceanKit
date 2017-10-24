scaleFactor = 1;
LoadFigureDefaults;

[rho, N2, zIn] = InternalModes.StratificationProfileWithName('exponential');
z = linspace(min(zIn),max(zIn),5000);
N0 = sqrt(max(N2(z)));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spectrum
%

omega = linspace(-N0,N0,200);
zOut = [0; -1300];

if isempty(GM)
    GM = GarrettMunkSpectrum(rho,zIn,latitude);
end
S_wkb = GM.HorizontalVelocitySpectrumAtFrequencies(zOut,omega, 'approximation', 'wkb');
S_gm = GM.HorizontalVelocitySpectrumAtFrequencies(zOut,omega, 'approximation', 'gm');
S = GM.HorizontalVelocitySpectrumAtFrequencies(zOut,omega);

S_wkb( S_wkb <1e-6 ) = 1e-6;
S_gm( S_gm <1e-6 ) = 1e-6;
S( S<1e-6 ) = 1e-6;

%%%%%%

xAxisMax = N0;

ticks = linspace(0,xAxisMax,5);
ticks = [-flip(ticks), ticks(2:end)];

labels = cell(length(ticks),1);
for i=1:length(ticks)
    if round(ticks(i)/f0) == 0.0
        labels{i} = '0';
    elseif round(ticks(i)/f0) == 1
        labels{i} = 'f_0';
    elseif round(ticks(i)/f0) == -1
        labels{i} = '-f_0';
    elseif ticks(i) == N0
        labels{i} = 'N_0';
    elseif ticks(i) == -N0
        labels{i} = '-N_0';
    else
        labels{i} = sprintf('%df_0',round(ticks(i)/f0));
    end
end

legendLabels = cell(length(zOut),1);
for i=1:length(zOut)
   if round(zOut(i)) == 0
       legendLabels{i} = 'surface';
   else
      legendLabels{i} = sprintf('%d m',round(zOut(i))); 
   end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Figure!
%
FigureSize = [50 50 figure_width_2col+8 175*scaleFactor];

fig1 = figure('Units', 'points', 'Position', FigureSize);
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');
fig1.PaperUnits = 'points';
fig1.PaperPosition = FigureSize;
fig1.PaperSize = [FigureSize(3) FigureSize(4)];

% subplot(1,2,1)
plot(omega,S, 'LineWidth', 1.0*scaleFactor), ylog, hold on
ax = gca;
ax.ColorOrderIndex = 1;
plot(omega,S_wkb,'LineStyle', '--', 'LineWidth', 1.0*scaleFactor)

ax.ColorOrderIndex = 1;
plot(omega,S_gm,'LineStyle', ':', 'LineWidth', 1.0*scaleFactor)

% ylim([1e0 1.1*max(max(S))])
ylim([1e-3 1e2])
ylabel('power ($m^2/s$)','Interpreter','LaTex', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
leg = legend(legendLabels);
% leg.Position(1) = 0.27;
% leg.Position(2) = 0.17;
set( gca, 'FontSize', figure_axis_tick_size);
xticks(ticks)
xticklabels(labels)
% set(gca, 'YTick', 1000*(-5:1:0));
% title('Horizontal velocity spectrum','Interpreter','LaTex', 'FontSize', figure_axis_label_size, 'FontName', figure_font);
