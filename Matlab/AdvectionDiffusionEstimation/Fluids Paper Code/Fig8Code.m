clear all
NX = 9; % number of drifters in X
dt = 1800; % time increment in seconds
NSim = 100; %Number of simulations
kappa_sm=.1;
%% Strain only
delta = 0;  
theta = 30;
zeta = 0; 
XX = 0.5:0.1:12;
rng(123)
    jjj = 1;
for NT = 1:300 %300; % number of time increments for simulation
    jj = 1;
for xx = XX
    [jjj xx]
    sigma = (1/xx)/(60*60*24);
    sigma_s = sigma.*sind(2*theta);
    sigma_n = sigma.*cosd(2*theta);
for j=1:NSim
    load("smoothedGriddedRho1Drifters.mat")
    x = x(1,1:NX); y = y(1,1:NX);
    for ii = 1:NT
        x(ii+1,:) = x(ii,:)+0.5*(sigma_n+delta).*x(ii,:)*dt+0.5*(sigma_s-zeta).*y(ii,:)*dt+sqrt(2*dt)*sqrt(kappa_sm)*randn(1,NX);
        y(ii+1,:) = y(ii,:)+0.5*(sigma_s+zeta).*x(ii,:)*dt+0.5*(delta-sigma_n).*y(ii,:)*dt+sqrt(2*dt)*sqrt(kappa_sm)*randn(1,NX);
    end
    xc = x-mean(x,2); yc = y-mean(y,2);
    u = diff(xc)/dt; v = diff(yc)/dt;
    xc = xc(1:size(x,1)-1,:); yc = yc(1:size(x,1)-1,:);
    [divest(j), vortest(j), nstrainest(j), sstrainest(j)] = DivVortStrainEst(xc,yc,u,v,dt,1,1,0,0);
    strainest(j)=sqrt(nstrainest(j).^2 + sstrainest(j).^2);
end
SEM(jjj,jj)=std(strainest)/sigma;
jj = jj+1;
end
jjj = jjj+1;
end
%%
XX = 0.5:0.1:12;
fig8a=figure; p=pcolor(XX,(2:300)/48,min(SEM,0.9));
xlabel('1/\sigma (days)'); ylabel('window length (days)');
colorbar
set(p, 'EdgeColor', 'none');
hold on; plot(XX,(sum(SEM>0.5)+1)/48,'r','linewidth',2)
ylim([2/48 6]);
exportfig(fig8a, 'figures/bootstrap.eps', 'width', 16, 'color', 'cmyk','Fontmode','fixed','FontSize', 14);
save('SEM1.mat','SEM');