sigma = 1e-5;
theta = 30*pi/180;
zeta = 0;
dt = 1800;
T = 5*86400;
kappa = 0.1;

x0 = 1000;
y0 = 1000;

% n = 2048;
% t = linspace(0,T,n).';
t = (0:dt:(T-dt)).';
n = length(t);

model = LinearVelocityField(sigma, theta, zeta);

% Without diffusivity
[x,y] = model.ParticlePath(x0,y0,t,0,0,0);
dt = t(2)-t(1);
v = diff(x)/dt + sqrt(-1)*diff(y)/dt;

% With diffusivity
integrator = AdvectionDiffusionIntegrator(model,kappa);
[t_s,x_s,y_s] = integrator.particleTrajectories(x0,y0,T-dt,dt);
v = diff(x_s)/dt + sqrt(-1)*diff(y_s)/dt;

[fn, sn] = powspec( dt, v );

r = sqrt(x0*x0+y0*y0);
alpha = atan2(y0,x0);

if sigma > abs(zeta)
    s = sqrt(sigma*sigma-zeta*zeta);
    A = sigma + zeta*sin(2*(theta-alpha));
    B = s*cos(2*(theta-alpha));
    C = zeta + sigma*sin(2*(theta-alpha));
    
    MSD = @(r,alpha) (2*r*r/(s*s*s*T))*sinh(s*T/2)*(sigma*A*cosh(s*T/2) + sigma*B*sinh(s*T/2)) - (r*r/(s*s))*zeta*C;
    S = @(f) 2*(r*r/(T)) * sinh(s*T/4).^2 * ( (sigma*A*cosh(s*T/2) + sigma*B*sinh(s*T/2) - zeta*C)./(s*s/4 + (2*pi*f).^2) + s*s*C*(zeta/2)./(s*s/4 + (2*pi*f).^2).^2 );
    
    Spp = @(f) (r*r/(T)) * sinh(s*T/4).^2 * ( (sigma*A*cosh(s*T/2) + sigma*B*sinh(s*T/2) - zeta*C)./(s*s/4 + (2*pi*f).^2) + s*s*C*((2*pi*f)+zeta/2)./(s*s/4 + (2*pi*f).^2).^2 );
    Snn = @(f) (r*r/(T)) * sinh(s*T/4).^2 * ( (sigma*A*cosh(s*T/2) + sigma*B*sinh(s*T/2) - zeta*C)./(s*s/4 + (2*pi*f).^2) + s*s*C*(-(2*pi*f)+zeta/2)./(s*s/4 + (2*pi*f).^2).^2 );
    Rc = @(f) -s*s*C*(2*pi*f) ./ ( (sigma*A*cosh(s*T/2) + sigma*B*sinh(s*T/2) - zeta*C).*(s*s/4 + (2*pi*f).^2) + s*s*zeta*C);
    TV = (r*r/(2*s*T))*sinh(s*T/2)*sigma*(A*cosh(s*T/2)+B*sinh(s*T/2)) + r*r*zeta*C/4;
    mean(abs(v).^2)
else
    s = sqrt(zeta*zeta-sigma*sigma);
    
    A = (1/2)*(zeta+s)*(zeta+sigma*sin(2*(theta-alpha)))*(sin(s*T/4)^2);
    B = (1/2)*(zeta-s)*(zeta+sigma*sin(2*(theta-alpha)))*(sin(s*T/4)^2);
    C = (sigma*(sigma+zeta*sin(2*(theta-alpha)))*cos(s*T/2) - sigma*s*cos(2*(theta-alpha))*sin(s*T/2))*(sin(s*T/4)^2);
    S = @(f) (r*r/T)*( A./((2*pi*f)-s/2).^2 + B./((2*pi*f)+s/2).^2 - C./(((2*pi*f)-s/2).*((2*pi*f)+s/2)) );
    
    % The Fourier transform of w = u+iv:
    Omega_p = @(f) (2*pi*f+s/2)*T/2;
    Omega_m = @(f) (2*pi*f-s/2)*T/2;
    w = @(f) (r/2)*( sigma*exp(2*sqrt(-1)*alpha) + sqrt(-1)*(zeta+s))*exp(-sqrt(-1)*Omega_m(f)) .* sin( Omega_m(f) )./(2*pi*f - s/2)  + (r/2)*( sigma*exp(2*sqrt(-1)*alpha) + sqrt(-1)*(zeta-s))*exp(-sqrt(-1)*Omega_p(f)) .* sin( Omega_p(f) )./(2*pi*f + s/2) ;
    
    A = sigma + zeta*sin(2*(theta-alpha));
    B = s*cos(2*(theta-alpha));
    C = zeta + sigma*sin(2*(theta-alpha));
    
    S = @(f) (r*r/T) * sin(s*T/4).^2 * ( -(sigma*A*cos(s*T/2) - sigma*B*sin(s*T/2) - zeta*C)./(-s*s/4 + (2*pi*f).^2) + s*s*C*((2*pi*f)+zeta/2)./(-s*s/4 + (2*pi*f).^2).^2 );
end

% MSD(sqrt(x0*x0+y0*y0),atan2(y0,x0))
% mean(x.^2 + y.^2)

snn = flip(sn(1:(floor(n/2)-1)));
f_snn = abs(flip(fn(1:(floor(n/2)-1))));
spp = sn((floor(n/2)+1):end);
f_spp = fn((floor(n/2)+1):end);

frequencyScale = 86400;

% We are plotting the *one-sided* spectrum. The theoretical value for S has
% been defined this way already, but this means we need double the
% theoretical white noise spectrum and add together spp and snn for the
% particle.

scaleFactor = 1;
LoadFigureDefaults

figure('Units', 'points', 'Position', [50 50 figure_width_1col 175*scaleFactor])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

plot( frequencyScale*f_spp, S(f_spp), 'LineWidth', 2 ), ylog, xlog
hold on, plot(frequencyScale*f_spp, 2*4*kappa*ones(size(f_spp)), 'LineWidth', 2)
plot(frequencyScale*f_spp,snn+spp, 'Color', 0.3*[1 1 1], 'LineWidth', 1.5)
xlim([frequencyScale*f_spp(1) frequencyScale*f_spp(end)])
ylim([1e-2 6e1])
xlabel('frequency (cycles/day)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
ylabel('power (m^2/s)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
p1 = gca;
p1.FontSize = figure_axis_tick_size;

print('SpectrumStrainOnlyJJE.eps','-depsc2')
