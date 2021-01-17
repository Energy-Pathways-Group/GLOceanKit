clear all
NX = 9; % number of drifters in X
NT = 49; %300; % number of time increments for simulation
dt = 1800; % time increment in seconds
NSim = 100; %Number of simulations
sigma=7e-6; theta=30;
PP = [0 0 0 0; sigma*cosd(2*theta) sigma*sind(2*theta) 0 0; sigma*cosd(2*theta) sigma*sind(2*theta) 6e-6 0; sigma*cosd(2*theta) sigma*sind(2*theta) 8e-6 0];
kappa_sm=.1;
Models = [1 1 0 0; 1 0 0 0; 0 1 0 0; 0 0 0 0];
T=dt*(NT-1);
f = -1/(2*dt):1/T:( 1/(2*dt) -1/T );
omega = 2*pi.*f;
for jj = 1:4
    rng(123)
    for kk = 1:NSim
    load("smoothedGriddedRho1Drifters.mat")
    x0 = x(1,1:NX); y0 = y(1,1:NX);
    x=zeros(NT,NX); y=zeros(NT,NX);
    x(1,:)=x0; y(1,:)=y0;  
    for ii = 1:NT
        x(ii+1,:) = x(ii,:)+0.5*(PP(jj,1)+PP(jj,4))*x(ii,:)*dt+0.5*(PP(jj,2)-PP(jj,3))*y(ii,:)*dt+sqrt(2*dt)*sqrt(kappa_sm)*randn(1,NX);
        y(ii+1,:) = y(ii,:)+0.5*(PP(jj,2)+PP(jj,3))*x(ii,:)*dt+0.5*(PP(jj,4)-PP(jj,1))*y(ii,:)*dt+sqrt(2*dt)*sqrt(kappa_sm)*randn(1,NX);
    end
    xc = x-mean(x,2); yc = y-mean(y,2);
    u = diff(xc)/dt; v = diff(yc)/dt;
    xc = xc(1:length(x)-1,:); yc = yc(1:length(y)-1,:);
    for ii = 1:size(Models,1)
        [delta, zeta, sigman, sigmas, u1, v1, u_sm, v_sm] = DivVortStrainEst( xc, yc, u, v, dt, Models(ii,1), Models(ii,2), Models(ii,3), Models(ii,4));
        FVU(jj,ii,kk) = mean(mean(u_sm.^2+v_sm.^2,2))./mean(mean(u.^2+v.^2,2));
        FDU(jj,ii,kk) = mean(abs(sum(u_sm+1i*v_sm,1)).^2)/mean(abs(sum(u+1i*v,1)).^2);
    end
        u_sm = u - 0.5*(PP(jj,1)+PP(jj,4))*xc - 0.5*(PP(jj,2)-PP(jj,3))*yc;
        v_sm = v - 0.5*(PP(jj,2)+PP(jj,3))*xc - 0.5*(PP(jj,4)-PP(jj,1))*yc;
        FVU(jj,ii+1,kk) = mean(mean(u_sm.^2+v_sm.^2,2))./mean(mean(u.^2+v.^2,2));
        FDU(jj,ii+1,kk) = mean(abs(sum(u_sm+1i*v_sm,1)).^2)/mean(abs(sum(u+1i*v,1)).^2);
    end
    if sum(PP(jj,:))==0
    FDUt(jj)=1; FVUt(jj)=1;
    else 
    zeta = PP(jj,3);
    thetar = theta*(pi/180);
    alpha = atan(yc(1,:)./xc(1,:)); %calc per particle
    r = sqrt(xc(1,:).^2+yc(1,:).^2); %calc per particle
    A = sigma + zeta*sin(2*(thetar-alpha));
    B = sigma*cos(2*(thetar-alpha));
    C = zeta + sigma*sin(2*(thetar-alpha));
        if zeta<sigma
            s = sqrt(sigma^2-zeta^2);
            S_fun = @(omega) r.^2./T.*sinh(s.*T/4).^2.* ( (sigma.*A.*cosh(s.*T./2)+sigma.*B.*sinh(s.*T./2)-zeta.*C)./(s.^2./4+omega.^2) + (s.^2.*C.*(omega+zeta/2))./(s.^2./4+omega.^2).^2 ) + (8/9)*4*kappa_sm;
        else
            s = sqrt(zeta^2-sigma^2);
            S_fun = @(omega) r.^2./T.*sin(s.*T/4).^2.* ( (-sigma.*A.*cos(s.*T./2)+sigma.*B.*sin(s.*T./2)+zeta.*C)./(-s.^2./4+omega.^2) + (s.^2.*C.*(omega+zeta/2))./(-s.^2./4+omega.^2).^2 ) + (8/9)*4*kappa_sm;
        end
    FDUt(jj) = 9*(8/9)*kappa_sm/sum(0.25*S_fun(0))
    FVUt(jj) = 4*(8/9)*kappa_sm/mean(mean(S_fun(omega')))
    end
end
%%
boxfig = figure;
for jj = 1:4
subplot(4,2,2*jj-1)
boxplot([squeeze(FVU(jj,1,:)), squeeze(FVU(jj,2,:)), squeeze(FVU(jj,3,:)), squeeze(FVU(jj,4,:)), squeeze(FVU(jj,5,:))], 'Labels',{'\{\sigma \theta\}', '\{\sigma \theta, \zeta\}',  '\{\sigma \theta, \delta\}', '\{\sigma \theta, \zeta, \delta\}','true'},'OutlierSize',1);
h = findobj(gca, 'type', 'text');
set(gca, 'TickLabelInterpreter', 'tex');
if jj == 1
    ylabel('\{\}')
    title('FVU')
    set(gca,'xtick',[]) 
elseif jj == 2
    ylabel('\{\sigma, \theta\}')
    set(gca,'xtick',[]) 
elseif jj == 3
    ylabel('\{\sigma, \theta > \zeta\}')
    set(gca,'xtick',[]) 
else
    ylabel('\{\sigma, \theta < \zeta\}')
end
ylim([0.85 1.01])
yline(FVUt(jj), 'r')
end
for jj = 1:4
subplot(4,2,2*jj)
boxplot([squeeze(FDU(jj,1,:)), squeeze(FDU(jj,2,:)), squeeze(FDU(jj,3,:)), squeeze(FDU(jj,4,:)), squeeze(FDU(jj,5,:))], 'Labels',{'\{\sigma \theta\}', '\{\sigma \theta, \zeta\}',  '\{\sigma \theta, \delta\}', '\{\sigma \theta, \zeta, \delta\}','true'},'OutlierSize',1);
h = findobj(gca, 'type', 'text');
set(gca, 'TickLabelInterpreter', 'tex');
set(gca, 'YAxisLocation', 'right')
ylim([0.00 1.05])
yline(FDUt(jj), 'r')
if jj == 1
title('FDU')
end
if jj < 4
    set(gca,'xtick',[]) 
end
end
ha=get(gcf,'children');
set(ha(1),'position',[.5 .1 .39 .19])
set(ha(2),'position',[.5 .3 .39 .19])
set(ha(3),'position',[.5 .5 .39 .19])
set(ha(4),'position',[.5 0.7 .39 .19])
set(ha(5),'position',[.1 .1 .39 .19])
set(ha(6),'position',[.1 .3 .39 .19])
set(ha(7),'position',[.1 .5 .39 .19])
set(ha(8),'position',[.1 0.7 .39 .19])
exportfig(boxfig, 'figures/FVUFDU.eps', 'width', 32, 'color', 'cmyk','Fontmode','fixed','FontSize', 14);