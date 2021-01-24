tic
clear all
NT = 48; %300; % number of time increments for simulation
dt = 1800; % time increment in seconds
NSim = 100; %Number of simulations
kappa_sm=0.1;
%% Strain only
rng(123)
delta = 0;  
sigma = 7e-6; 
theta = 30;
zeta = 6e-6; 
sigma_s = sigma.*sind(2*theta);
sigma_n = sigma.*cosd(2*theta);
    load("smoothedGriddedRho1Drifters.mat")
    xl = x(1,1:9); yl = y(1,1:9);
    xlc = mean(xl); ylc = mean(yl);
    cmd = mean(sqrt((xl-xlc).^2+(yl-ylc).^2));
KK = 2:48;
profile on
for NX = KK
    NX
for j=1:NSim
    x = xlc+cmd*randn(1,NX); y = ylc+cmd*randn(1,NX);    
    for ii = 1:NT
        x(ii+1,:) = x(ii,:)+0.5*(sigma_n+delta).*x(ii,:)*dt+0.5*(sigma_s-zeta).*y(ii,:)*dt+sqrt(2*dt)*sqrt(kappa_sm)*randn(1,NX);
        y(ii+1,:) = y(ii,:)+0.5*(sigma_s+zeta).*x(ii,:)*dt+0.5*(delta-sigma_n).*y(ii,:)*dt+sqrt(2*dt)*sqrt(kappa_sm)*randn(1,NX);
    end
    xc = x-mean(x,2); yc = y-mean(y,2);
    u = diff(xc)/dt; v = diff(yc)/dt;
    xc = xc(1:length(x)-1,:); yc = yc(1:length(y)-1,:);
    if 1 == 1
        [divest(j), vortest(j), nstrainest(j), sstrainest(j)] = DivVortStrainEst(xc,yc,u,v,dt,1,0,0,0);
    else
        parameterEstimates = EstimateLinearVelocityFieldParameters( xc, yc, dt*(1:NT)', [ModelParameter.strain,ModelParameter.vorticity] );
        divest(j) = parameterEstimates.delta;
        vortest(j) = parameterEstimates.zeta;
        nstrainest(j) = parameterEstimates.sigma_n;
        sstrainest(j) = parameterEstimates.sigma_s;
    end
    strainest(j)=sqrt(nstrainest(j).^2 + sstrainest(j).^2);
    thetaest(j)=atan2d(sstrainest(j),nstrainest(j))/2;
end
%%
STD(NX-1,:) = [std(strainest)/sigma std(thetaest)/theta std(vortest)/zeta];
end
profile viewer
%%
figure; plot(KK,STD)
figa0=figure; loglog(KK,STD,'linewidth',2); hold on; loglog(5:45, .3./sqrt(5:45),'linewidth',2)
legend('strain rate (\sigma)','strain angle (\theta)','vorticity (\zeta)','K^{-1/2}')
xlabel('number of drifters (K)'); ylabel('relative standard error')
xlim([min(KK) max(KK)])
exportfig(figa0, 'figures/append0.eps', 'width', 16, 'color', 'cmyk','Fontmode','fixed','FontSize', 14);
toc