scaleFactor = 1;
LoadFigureDefaults;

latitude = 33;
f0 = 2 * 7.2921E-5 * sin( latitude*pi/180 );
L = 5000;

j_star = 3;
L_gm = 1.3e3; % thermocline exponential scale, meters
invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_gm = 6.3e-5; % non-dimensional energy parameter
g = 9.81;
rho0 = 1025;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Constant stratification structure function
%

N0 = invT_gm/4; % 'average' GM bouyancy frequency, sqrt( (1-exp(-L/L_gm))/(exp(L/L_gm-1) -1))
rho = @(z) -(N0*N0*rho0/g)*z + rho0;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Exponential stratification structure function
%
N0 = invT_gm;
rho = @(z) rho0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*z/L_gm)));


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Realistic stratification structure function
%
% delta_p = 75;
% z_s = 200;
% L_s = 500;
% L_d = 2500; %L_d = 1300;
% z_p = -500; %z_p = -1;
% D = 5000;
% 
% N0 = 0.5e-3; %N0 = 5.2e-3;
% Np = 10e-3;
% Nd = 3e-4; %Nd = 1.1e-4;
% 
% A = Np*Np - Nd*Nd*exp(2*(D+z_p)/L_d);
% B = (Nd*Nd*exp(2*(D+z_p)/L_d) - N0*N0)/(exp(-2*(z_p-z_s)/L_s) - exp(2*z_s/L_s));
% C = N0*N0-B*exp(2*z_s/L_s);
% E = Nd*Nd*exp(2*(D+z_p)/L_d);
% 
% 
% rho_surface = @(z) rho0*(1 - (L_s*B/(2*g)) * (exp(2*z_s/L_s) - exp(-2*(z-z_s)/L_s)) - C*z/g);
% rho_deep = @(z) rho_surface(z_p) + (rho0*L_d*E/(2*g))*(1-exp(2*(z-z_p)/L_d));
% rho_p = @(z) -(A*rho0*delta_p/g)*(tanh( (z-z_p)/delta_p) - 1);
% 
% rho = @(z) (z>=z_p).*rho_surface(z) + (z<z_p).*rho_deep(z) + rho_p(z);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Spectrum
%

omega = linspace(-N0,N0,200);
% zOut = [0; -L_gm/2; -L_gm];
zOut = [0; -L_gm];

GM = GarrettMunkSpectrum(rho,[-L 0],latitude);
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
