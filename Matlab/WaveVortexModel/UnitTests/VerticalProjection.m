% Playing with convergence properties of the vertical mode projection.
%
% 2021-11-14
% If you anti-alias, and therefore truncate the modes, then in physical
% space oscillations become apparent. Essentially you're using the missing
% variance that would have been captured by the next few modes. The silly
% ah-ha that I had, was that of course, when you're using the quadrature
% grid, and haven't truncated, you won't *see* the oscillations, but
% they are, of course, still there. So really the only solution is to go to
% higher vertical resolution. That's it.

Lx = 750e3;
Ly = 750e3;

N = 16;
Nx = N;
Ny = N;
nModes = 75;

latitude = 31;
f0 = 2 * 7.2921E-5 * sin( latitude*pi/180 );

N0 = 12*2*pi/3600; % reference buoyancy frequency, radians/seconds
Nmin = N0*3e-2;
rho0 = 1025; g = 9.81;
L_gm = 145; % thermocline exponential scale, meters
rhoFunction = @(z) rho0*(1 + L_gm*N0*N0/(2*g)*(1 - exp(2*z/L_gm)) - (Nmin*Nmin/g)*z);
N2Function = @(z) N0*N0*exp(2*z/L_gm) + Nmin*Nmin*ones(size(z));
dLnN2Function = @(z) 2*ones(size(z))/L_gm;

zIn = [-2000 0];

A = -1.32e-5; % s^{-1}
alpha = 8e-10; % m^{-2}
beta = 8.2e-6; % m^{-2}
etaFunction = @(z) -A*f0*(beta/alpha)*(z./N2Function(z)).*exp(-beta*z.*z) - A*A*(beta/alpha)*(z./N2Function(z)).*exp(-2*beta*z.*z);


% %%%%%%%%%%%%%%%%%%%%%%%%%
% % small
% wvmSmall = WaveVortexModelHydrostatic([Lx, Ly, max(zIn)-min(zIn)], [Nx, Ny, nModes/2], latitude, rhoFunction,'N2func', N2Function, 'dLnN2func',dLnN2Function);
% eta = etaFunction(wvmSmall.z);
% 
% Q = squeeze(wvmSmall.Q);
% sGinv = wvmSmall.QGinv.*shiftdim(Q,-1);
% sG = wvmSmall.QG./Q;
% 
% etahat_s = sG*eta;

%%%%%%%%%%%%%%%%%%%%%%%%%
% full

for nModes = 100 %20:10:100
wvm = WaveVortexModelHydrostatic([Lx, Ly, max(zIn)-min(zIn)], [Nx, Ny, nModes], latitude, rhoFunction,'N2func', N2Function, 'dLnN2func',dLnN2Function);



% operators
Q = squeeze(wvm.Q);
Ginv = wvm.QGinv.*shiftdim(Q,-1);
G = wvm.QG./Q;

% filter

minRMS = 1e9;
% mint, minp, minalphaf;
for t=(1:10)/10
    for p=2:2:8
        for alphaf=(1:8)/1
%             t=7/8;p=8;alphaf=0.5;
            j_max = max(wvm.j);
            j_max = 2*j_max/3;
            dj = wvm.j(2)-wvm.j(1);
            j_cutoff = dj*(j_max/dj)^t;
            alpha = alphaf*log(10); % orders of magnitude to drop over the filter
            % p = 4;
            Q = @(x,x_max,x_cutoff) exp( - alpha*((abs(x)-x_cutoff)./(x_max-x_cutoff)).^p );
            Qj = Q(wvm.j,j_max,j_cutoff);
            Qj(isnan(Qj)) = 1;
            Qj(abs(wvm.j) < j_cutoff) = 1;
            Qj(abs(wvm.j) > j_max) = 0;

            eta = etaFunction(wvm.z);
            etahat = G*eta;

            etaback = Ginv*(etahat.*Qj);
            
%             rms = mean( (eta-etaback).^2);
            rms = mean(sqrt(diff(eta-etaback).^2));
%             fprintf('t=%.2f,p=%d,alpha=%.1f,rms=%.2f\n',t,p,alphaf,rms)
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
            j_cutoff = dj*(j_max/dj)^mint;
            alpha = minalphaf*log(10); % orders of magnitude to drop over the filter
            % p = 4;
            Q = @(x,x_max,x_cutoff) exp( - alpha*((abs(x)-x_cutoff)./(x_max-x_cutoff)).^minp );
            Qj = Q(wvm.j,j_max,j_cutoff);
            Qj(isnan(Qj)) = 1;
            Qj(abs(wvm.j) < j_cutoff) = 1;
            Qj(abs(wvm.j) > j_max) = 0;

            eta = etaFunction(wvm.z);
            etahat = G*eta;

            etaback = Ginv*(etahat.*Qj);

fprintf('t=%.2f (j_c=%.2f of %.2f),p=%d,alpha=%.4f,rms=%.2f\n',mint,j_cutoff, j_max, minp,minalphaf,minRMS)
figure, plot(eta,wvm.z), hold on, plot(etaback,wvm.z)
end