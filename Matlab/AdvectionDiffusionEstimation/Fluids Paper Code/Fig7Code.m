clear all
NX = 9; % number of drifters in X
NT = 300; %300; % number of time increments for simulation
dt = 1800; % time increment in seconds
NSim = 100; %Number of simulations
sigma = 1e-5:(1e-6-1e-5)/NT:(1e-6-(1e-6-1e-5)/NT); 
theta=30;
sigma_s = sigma*sind(2*theta);
sigma_n = sigma*cosd(2*theta);
kappa_sm=.5;
W = [12 48 144];
for ww = 1:3
rng(456)
for kk = 1:NSim
    load("smoothedGriddedRho1Drifters.mat")
    x0 = x(1,1:NX); y0 = y(1,1:NX);
    x=zeros(NT,NX); y=zeros(NT,NX);
    x(1,:)=x0; y(1,:)=y0;  
    for ii = 1:NT
        x(ii+1,:) = x(ii,:)+0.5*(sigma_n(ii))*x(ii,:)*dt+0.5*(sigma_s(ii))*y(ii,:)*dt+sqrt(2*dt)*sqrt(kappa_sm)*randn(1,NX);
        y(ii+1,:) = y(ii,:)+0.5*(sigma_s(ii))*x(ii,:)*dt+0.5*(-sigma_n(ii))*y(ii,:)*dt+sqrt(2*dt)*sqrt(kappa_sm)*randn(1,NX);
    end
    xc = x-mean(x,2); yc = y-mean(y,2);
    u = diff(xc)/dt; v = diff(yc)/dt;
    xc = xc(1:length(x)-1,:); yc = yc(1:length(y)-1,:);
    for jj = W(ww)/2+1:NT-1-W(ww)/2
        [a , b , sigman(ww,kk,jj), sigmas(ww,kk,jj), u1, v1, u_sm, v_sm] = DivVortStrainEst( xc(jj-W(ww)/2:jj+W(ww)/2,:), yc(jj-W(ww)/2:jj+W(ww)/2,:), u(jj-W(ww)/2:jj+W(ww)/2,:), v(jj-W(ww)/2:jj+W(ww)/2,:), dt, 1 ,1,0,0);
    end
end
end
%%
sigmaE = sqrt(sigman.^2+sigmas.^2);
StrainVar=figure; set(gcf,'Position',[100 100 2100 500])
subplot(1,2,1)
plot(dt*(1:NT)/86400,sigma,'k','linewidth',3); xlim([0 6.25]);
newcolors = [0.83 0.14 0.14; 1.00 0.54 0.00; 0.47 0.25 0.80]; colororder(newcolors)
for ii = 1:3
hold on; plot(t(W(ii)/2+1:NT-1-W(ii)/2)/86400,squeeze(sigmaE(ii,1,W(ii)/2+1:NT-1-W(ii)/2)),'linewidth',1.5)
end
legend('truth','W=6 hours','W = 1 day','W = 3 days');
xlabel('days'); ylabel('\sigma'); 
subplot(1,2,2); plot(dt*(1:NT)/86400,sigma/2,'k','linewidth',3); xlim([0 6.25]);
for ii = 1:3
hold on; plot(t(W(ii)/2+1:NT-1-W(ii)/2)/86400,(std(squeeze(sigmaE(ii,:,W(ii)/2+1:NT-1-W(ii)/2)))),'linewidth',2)
end
legend('\sigma/2','W=6 hours','W = 1 day','W = 3 days');
xlabel('days');  ylabel('standard error');
exportfig(StrainVar, 'figures/StrainVar.eps', 'width', 32, 'color', 'cmyk','Fontmode','fixed','FontSize', 14);