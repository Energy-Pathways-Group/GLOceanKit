% Playing with convergence properties of the vertical mode projection.
%
% 2021-11-14
% If you anti-alias, and therefore truncate the modes, then in physical
% space oscillations become apparent. Essentially you're losing the missing
% variance that would have been captured by the next few modes. The silly
% ah-ha that I had, was that of course, when you're using the quadrature
% grid, and haven't truncated, you won't *see* the oscillations, but
% they are, of course, still there. So really the only solution is to go to
% higher vertical resolution.
%
% 2021-11-15
% The projection of *eta* results is significantly larger oscillations in
% *rho*. The reason is that they differ by a factor of N^2 (which is orders
% of magnitude larger at the surface than at depth). So to keep a smooth
% density structure we do, in fact, want to filter the eta-modes.
%
% Looks like the the number of modes that need to be damped
% (j_max-j_cutoff) is the one thing that really needs to be varied somewhat
% manually.

Lx = 750e3;
Ly = 750e3;

N = 16;
Nx = N;
Ny = N;

latitude = 31;
f = 2 * 7.2921E-5 * sin( latitude*pi/180 );

% Cyprus Eddy example stratification profile
N0 = 12*2*pi/3600; % reference buoyancy frequency, radians/seconds
rho0 = 1025; g = 9.81;
L_gm = 145; % thermocline exponential scale, meters
L_const = 3*L_gm; % depth below which stratification stays constant--try 8*L_gm for badness, 3*L_gm for goodness
Nmin = sqrt(N0*N0*exp(-2*L_const/L_gm));
rhoFunction = @(z) rho0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*z/L_gm)) - (Nmin*Nmin/g)*z);
N2Function = @(z) N0*N0*exp(2*z/L_gm) + Nmin*Nmin*ones(size(z));
dLnN2Function = @(z) 2*ones(size(z))/L_gm;

zIn = [-2000 0];

% Cyprus Eddy buoyancy anomaly.
A = -1.32e-5; % s^{-1}
alpha = 8e-10; % m^{-2}
beta = 8.2e-6; % m^{-2}
% etaFunction = @(z) -A*f*(beta/alpha)*(z./N2Function(z)).*exp(-beta*z.*z) - A*A*(beta/alpha)*(z./N2Function(z)).*exp(-2*beta*z.*z);
eta = @(x,y,z) -A*f*(beta/alpha)*(z./N2Function(z)).*exp(-alpha*((x).*(x)+(y).*(y))-beta*z.*z) - A*A*(beta/alpha)*(z./N2Function(z)).*exp(-2*alpha*((x).*(x)+(y).*(y))-2*beta*z.*z);
etaFunction= @(z) eta(-3e3,-3e3,z);

% Start by plotting the *density* profile (not eta)
z = linspace(min(zIn),max(zIn),1000)';
figure
sp1 = subplot(2,2,1);
plot((rho0/g)*N2Function(z).*etaFunction(z),z); hold on
title('density anomaly')
ylabel('raw, anti-alaised')
sp2 = subplot(2,2,3);
plot((rho0/g)*N2Function(z).*etaFunction(z),z); hold on
xlabel('kg m^{-3}')
ylabel('filtered')
sp3 = subplot(2,2,2);
plot(etaFunction(z),z); hold on
title('eta')
sp4 = subplot(2,2,4);
plot(etaFunction(z),z); hold on
xlabel('m')

error = @(rho,rhoback) mean(sqrt(diff(rho-rhoback).^2));
error = @(rho,rhoback) sqrt(mean( (rho-rhoback).^2));

% Now loop through models with an increasing number of modes,
nModes = (20:10:100)';
j_cutoff_optimal = zeros(size(nModes));
j_max_optimal = zeros(size(nModes));
for iMode=1:length(nModes)
    wvm = WaveVortexModelHydrostatic([Lx, Ly, max(zIn)-min(zIn)], [Nx, Ny, nModes(iMode)], latitude, rhoFunction,'N2func', N2Function, 'dLnN2func',dLnN2Function);
    wvm.shouldAntiAlias = 1;
    
    % Vertical projection operators
    Q = squeeze(wvm.Q);
    Ginv = wvm.QGinv.*shiftdim(Q,-1);
    G = wvm.QG./Q;

    j_max = 2*max(wvm.j)/3;
    j_max_optimal(iMode) = j_max;
    
    % Now first, let's look at the baseline
    eta = etaFunction(wvm.z);
    rho = (rho0/g)*N2Function(wvm.z).*eta;
    etahat = G*eta;
    
    Qj = ones(size(etahat));
    Qj(round(j_max):end)=0;
    etaback = Ginv*(etahat.*Qj);
    rhoback = (rho0/g)*N2Function(wvm.z).*etaback;
    subplot(sp1), plot(rhoback,wvm.z)
    subplot(sp3), plot(etaback,wvm.z)
    
    baselineError = error(rho,rhoback);

    % Now let's minimize an error function to find the "optimal" spectral
    % filter
    minRMS = 1e9;
    % mint, minp, minalphaf;
    for t=(1:100)/100
        for p=8%2:2:8
            for alphaf=1%(1:15)/5
                %             t=7/8;p=8;alphaf=0.5;
                j_max = max(wvm.j);
                j_max = 2*j_max/3;
                dj = wvm.j(2)-wvm.j(1);
                j_cutoff = t*j_max; %dj*(j_max/dj)^t;
                alpha = alphaf*log(10); % orders of magnitude to drop over the filter
                Q = @(x,x_max,x_cutoff) exp( - alpha*((abs(x)-x_cutoff)./(x_max-x_cutoff)).^p );
                Qj = Q(wvm.j,j_max,j_cutoff);
                Qj(isnan(Qj)) = 1;
                Qj(abs(wvm.j) < j_cutoff) = 1;
                Qj(abs(wvm.j) > j_max) = 0;

                eta = etaFunction(wvm.z);
                rho = (rho0/g)*N2Function(wvm.z).*eta;
                etahat = G*eta;

                etaback = Ginv*(etahat.*Qj);
                rhoback = (rho0/g)*N2Function(wvm.z).*etaback;

                rms = error(rho,rhoback);

                if rms < minRMS
                    minRMS = rms;
                    mint = t; minp = p; minalphaf = alphaf;
                end
            end
        end
    end

    j_max = max(wvm.j);
    j_max = 2*j_max/3;
    dj = wvm.j(2)-wvm.j(1);
    j_cutoff = mint*j_max; % dj*(j_max/dj)^mint;
    alpha = minalphaf*log(10); % orders of magnitude to drop over the filter
    Q = @(x,x_max,x_cutoff) exp( - alpha*((abs(x)-x_cutoff)./(x_max-x_cutoff)).^minp );
    Qj = Q(wvm.j,j_max,j_cutoff);
    Qj(isnan(Qj)) = 1;
    Qj(abs(wvm.j) < j_cutoff) = 1;
    Qj(abs(wvm.j) > j_max) = 0;

    eta = etaFunction(wvm.z);
    rho = (rho0/g)*N2Function(wvm.z).*eta;
    etahat = G*eta;

    etaback = Ginv*(etahat.*Qj);
    rhoback = (rho0/g)*N2Function(wvm.z).*etaback;
    subplot(sp2), plot(rhoback,wvm.z)
    subplot(sp4), plot(etaback,wvm.z)
    
    j_cutoff_optimal(iMode) = j_cutoff;

    % minRMS = mean(sqrt(diff(eta-etaback).^2));
    % minRMS = mean(sqrt(diff(rho-rhoback).^2));

    fprintf('t=%.2f (j_c=%.2f of %.2f),p=%d,alpha=%.4f, error_ratio=%.2g\n',mint,j_cutoff, j_max, minp,minalphaf,minRMS/baselineError);
end
packfig(2,2)