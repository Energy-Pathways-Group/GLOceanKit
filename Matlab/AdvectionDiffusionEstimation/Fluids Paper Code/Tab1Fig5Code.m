clear all
NX = 9; % number of drifters in X
NT = 48; %300; % number of time increments for simulation
dt = 1800; % time increment in seconds
NSim = 100; %Number of simulations
kappa_sm=0.1;
%%
B=100;
SE_B = @(x,B) sqrt( (1/(B-1))*sum((x-mean(x,2)).^2,2) );
%% Strain only
rng(123)
delta = 0;  
sigma = 7e-6; 
theta = 30;
zeta = 0; 
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
    xc = xc(1:length(x)-1,:); yc = yc(1:length(y)-1,:);
    [divest(j), vortest(j), nstrainest(j), sstrainest(j)] = DivVortStrainEst(xc,yc,u,v,dt,1,1,0,0);
    strainest(j)=sqrt(nstrainest(j).^2 + sstrainest(j).^2);
    thetaest(j)=atan2d(sstrainest(j),nstrainest(j))/2;
    for i=1:B %create bootstrap samples
        num = datasample([1:size(xc,2)],size(xc,2),2);
        xB(:,:,i) = xc(:,num);
        yB(:,:,i) = yc(:,num);
        uB(:,:,i) = u(:,num);
        vB(:,:,i) = v(:,num);
    end
    for b = 1:size(xB,3) %loop over bootstrap samples
        x1 = xB(:,:,b);
        y1 = yB(:,:,b);
        u1 = uB(:,:,b);
        v1 = vB(:,:,b);
        [divestB(j,b), vortestB(j,b), nstrainestB(j,b), sstrainestB(j,b)] = DivVortStrainEst(x1,y1,u1,v1,dt,1,1,0,0);
    end
    strainestB(j,:) = sqrt(nstrainestB(j,:).^2 + sstrainestB(j,:).^2);
    thetaestB(j,:) = atan2d(sstrainestB(j,:),nstrainestB(j,:))/2;
    SE_strain_B(j) = SE_B(strainestB(j,:),B);
    SE_theta_B(j) = SE_B(thetaestB(j,:),B);
end
%%
SEfun = @(x) sqrt(var(x));
[SEfun(strainest); mean(SE_strain_B); SEfun(SE_strain_B)]
[SEfun(thetaest); mean(SE_theta_B); SEfun(SE_theta_B)]
%% Strain dominated
rng(123)
delta = 0;  
sigma = 7e-6; 
theta = 30;
zeta = 6e-6;
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
    xc = xc(1:length(x)-1,:); yc = yc(1:length(y)-1,:);
    [divest(j), vortest(j), nstrainest(j), sstrainest(j)] = DivVortStrainEst(xc,yc,u,v,dt,1,0,0,0);
    strainest(j)=sqrt(nstrainest(j).^2 + sstrainest(j).^2);
    thetaest(j)=atan2d(sstrainest(j),nstrainest(j))/2;
    for i=1:B %create bootstrap samples
        num = datasample([1:size(xc,2)],size(xc,2),2);
        xB(:,:,i) = xc(:,num);
        yB(:,:,i) = yc(:,num);
        uB(:,:,i) = u(:,num);
        vB(:,:,i) = v(:,num);
    end
    for b = 1:size(xB,3) %loop over bootstrap samples
        x1 = xB(:,:,b);
        y1 = yB(:,:,b);
        u1 = uB(:,:,b);
        v1 = vB(:,:,b);
        [divestB(j,b), vortestB(j,b), nstrainestB(j,b), sstrainestB(j,b)] = DivVortStrainEst(x1,y1,u1,v1,dt,1,0,0,0);    
    end
    strainestB(j,:) = sqrt(nstrainestB(j,:).^2 + sstrainestB(j,:).^2);
    thetaestB(j,:) = atan2d(sstrainestB(j,:),nstrainestB(j,:))/2;
    SE_strain_B(j) = SE_B(strainestB(j,:),B);
    SE_theta_B(j) = SE_B(thetaestB(j,:),B);
    SE_vort_B(j) = SE_B(vortestB(j,:),B);
end
%%
SEfun = @(x) sqrt(var(x));
[SEfun(strainest); mean(SE_strain_B); SEfun(SE_strain_B)]
[SEfun(thetaest); mean(SE_theta_B); SEfun(SE_theta_B)]
[SEfun(vortest); mean(SE_vort_B); SEfun(SE_vort_B)]
%%
fig1 = figure; set(gcf,'Position',[100 100 2100 500])
subplot(1,3,1)
histogram(min(strainestB(:),2*sigma),50); l1 = xline(sigma,'r','linewidth',2); l2 = xline(mean(strainestB(:)),'b','linewidth',2);
xlim([0 2*sigma]); xlabel('\sigma')
subplot(1,3,2)
histogram(max(min(thetaestB(:),90),-30),50); l1 = xline(theta,'r','linewidth',2); l2 = xline(mean(thetaestB(:)),'b','linewidth',2);
xlim([-30 90]); xlabel('\theta')
subplot(1,3,3)
histogram(max(min(vortestB(:),3*zeta),-zeta),50); l1 = xline(zeta,'r','linewidth',2); l2 = xline(mean(vortestB(:)),'b','linewidth',2);
xlim([-zeta 3*zeta]); xlabel('\zeta')
exportfig(fig1, 'figures/bootstrap_histStrainDom.eps', 'width', 32, 'color', 'cmyk','Fontmode','fixed','FontSize', 14);