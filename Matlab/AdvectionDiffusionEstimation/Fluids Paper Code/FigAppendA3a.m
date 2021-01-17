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
for NT = 2:2:300 %300; % number of time increments for simulation
    jj = 1;
for xx = XX
    [jjj xx]
    sigma = (1/xx)/(60*60*24);
    sigma_s = sigma.*sind(2*theta);
    sigma_n = sigma.*cosd(2*theta);
for j=1:NSim
    load("smoothedGriddedRho1Drifters.mat")
    xl = x(1,1:NX); yl = y(1,1:NX);
    x(5) = mean(xl); y(5) = mean(yl);
    cmd = mean(sqrt((xl-x(5)).^2+(yl-y(5)).^2));
    spc = 9*cmd/20;
    x = (-4:4)*spc*cosd(30)+x(5); y = (-4:4)*spc*sind(30)+y(5);
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
save('SEM3.mat','SEM');

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
for NT = 2:2:300 %300; % number of time increments for simulation
    jj = 1;
for xx = XX
    [jjj xx]
    sigma = (1/xx)/(60*60*24);
    sigma_s = sigma.*sind(2*theta);
    sigma_n = sigma.*cosd(2*theta);
for j=1:NSim
    load("smoothedGriddedRho1Drifters.mat")
    xl = x(1,1:NX); yl = y(1,1:NX);
    x(5) = mean(xl); y(5) = mean(yl);
    cmd = mean(sqrt((xl-x(5)).^2+(yl-y(5)).^2));
    spc = 9*cmd/20;
    x = (-4:4)*spc*cosd(120)+x(5); y = (-4:4)*spc*sind(120)+y(5);
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
save('SEM4.mat','SEM');