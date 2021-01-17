clear all
load('smoothedGriddedRho1Drifters.mat')
ND = 9;
W = 48;
dt = t(2)-t(1);
f0 = 2*(7.2921e-5)*sin(lat0*pi/180);
x=x(:,1:ND)-mean(x(:,1:ND),2); y=y(:,1:ND)-mean(y(:,1:ND),2);
u = diff(x(:,1:ND))/dt; v = diff(y(:,1:ND))/dt;
x1 = x(1:length(x)-1,1:ND); y1 = y(1:length(y)-1,1:ND);
    for jj = W/2+1:length(t)-1-W/2
        [delta(jj), zeta(jj), sigman(jj), sigmas(jj)] = DivVortStrainEst( x1(jj-W/2:jj+W/2,:), y1(jj-W/2:jj+W/2,:), u(jj-W/2:jj+W/2,:), v(jj-W/2:jj+W/2,:), dt, 1,1,0,0);
    end
sigma = sqrt(sigman.^2+sigmas.^2)/f0;
theta = atan2d(sigmas,sigman)/2;
%%
B=100;
for b=1:B %create bootstrap samples
    num = datasample([1:size(x,2)],size(x,2),2);
    xB = x(:,num);
    yB = y(:,num);
    ub = diff(xB(:,1:ND))/dt; vb = diff(yB(:,1:ND))/dt;
    xb = xB(1:length(x)-1,1:ND); yb = yB(1:length(y)-1,1:ND);
    for jj = W/2+1:length(t)-1-W/2
        [deltab(b,jj), zetab(b,jj), sigmanb(b,jj), sigmasb(b,jj)] = DivVortStrainEst( xb(jj-W/2:jj+W/2,:), yb(jj-W/2:jj+W/2,:), ub(jj-W/2:jj+W/2,:), vb(jj-W/2:jj+W/2,:), dt, 1,1,0,0);
    end
end
sigmab = sqrt(sigmanb.^2+sigmasb.^2)/f0;
thetab = atan2d(sigmasb,sigmanb)/2;
%%
latmixroll = figure; 
subplot(2,2,1);
plot(t(W/2+1:length(t)-1-W/2)/86400,sigmab(:,W/2+1:length(t)-1-W/2),'color',[.5 .5 .5]);
hold on; plot(t/86400,0.0591*ones(1,length(t)),'r','linewidth',2); xlim([t(1)/86400 t(end)/86400])
hold on; plot(t(W/2+1:length(t)-1-W/2)/86400,sigma(W/2+1:length(t)-1-W/2),'b','linewidth',2)
xlabel('days'); ylabel('\sigma (f_0)');  title('Site 1 (strain model)')
%%
load('smoothedGriddedRho2Drifters.mat')
ND = 8;
W = 48;
dt = t(2)-t(1);
f0 = 2*(7.2921e-5)*sin(lat0*pi/180);
x=x(:,1:ND)-mean(x(:,1:ND),2); y=y(:,1:ND)-mean(y(:,1:ND),2);
u = diff(x(:,1:ND))/dt; v = diff(y(:,1:ND))/dt;
x1 = x(1:length(x)-1,1:ND); y1 = y(1:length(y)-1,1:ND);
Models = [1 1 0 0; 1 0 0 0];
for ii = 1:size(Models,1)
    for jj = W/2+1:length(t)-1-W/2
        [delta2(ii,jj), zeta2(ii,jj), sigman2(ii,jj), sigmas2(ii,jj), u1, v1, u_sm, v_sm] = DivVortStrainEst( x1(jj-W/2:jj+W/2,:), y1(jj-W/2:jj+W/2,:), u(jj-W/2:jj+W/2,:), v(jj-W/2:jj+W/2,:), dt, Models(ii,1), Models(ii,2), Models(ii,3), Models(ii,4));
    end
end
sigma = sqrt(sigman2.^2+sigmas2.^2)/f0;
theta = atan2d(sigmas2,sigman2)/2;
%%
B=100;
for ii = 1:size(Models,1)
for b=1:B %create bootstrap samples
    num = datasample([1:size(x,2)],size(x,2),2);
    xB = x(:,num);
    yB = y(:,num);
    ub = diff(xB(:,1:ND))/dt; vb = diff(yB(:,1:ND))/dt;
    xb = xB(1:length(x)-1,1:ND); yb = yB(1:length(y)-1,1:ND);
    for jj = W/2+1:length(t)-1-W/2
        [delta2b(ii,b,jj), zeta2b(ii,b,jj), sigman2b(ii,b,jj), sigmas2b(ii,b,jj)] = DivVortStrainEst( xb(jj-W/2:jj+W/2,:), yb(jj-W/2:jj+W/2,:), ub(jj-W/2:jj+W/2,:), vb(jj-W/2:jj+W/2,:), dt, Models(ii,1), Models(ii,2), Models(ii,3), Models(ii,4));
    end
end
end
sigmab = sqrt(sigman2b.^2+sigmas2b.^2)/f0;
thetab = atan2d(sigmas2b,sigman2b)/2;
%%
subplot(2,2,2);
plot(t(W/2+1:length(t)-1-W/2)/86400,reshape(sigmab(1,:,W/2+1:length(t)-1-W/2),B,255),'color',[.5 .5 .5]);
hold on; plot(t/86400,0.0131*ones(1,length(t)),'r','linewidth',2); xlim([t(1)/86400 t(end)/86400])
hold on; plot(t(W/2+1:length(t)-1-W/2)/86400,sigma(1,W/2+1:length(t)-1-W/2),'b','linewidth',2)
xlabel('days'); ylabel('\sigma (f_0)'); title('Site 2 (strain model)')
%%
subplot(2,2,3);
plot(t(W/2+1:length(t)-1-W/2)/86400,reshape(sigmab(2,:,W/2+1:length(t)-1-W/2),B,255),'color',[.5 .5 .5]);
hold on; plot(t/86400,0.0642*ones(1,length(t)),'r','linewidth',2); xlim([t(1)/86400 t(end)/86400])
hold on; plot(t(W/2+1:length(t)-1-W/2)/86400,sigma(2,W/2+1:length(t)-1-W/2),'b','linewidth',2)
xlabel('days'); ylabel('\sigma (f_0)');  title('Site 2 (strain-vorticity model)')
%%
subplot(2,2,4);
plot(t(W/2+1:length(t)-1-W/2)/86400,reshape(zeta2b(2,:,W/2+1:length(t)-1-W/2),B,255)/f0,'color',[.5 .5 .5]);
hold on; plot(t/86400,0.0637*ones(1,length(t)),'r','linewidth',2); xlim([t(1)/86400 t(end)/86400])
hold on; plot(t(W/2+1:length(t)-1-W/2)/86400,zeta2(2,W/2+1:length(t)-1-W/2)/f0,'b','linewidth',2)
xlabel('days'); ylabel('\zeta (f_0)');  title('Site 2 (strain-vorticity model)')
%%
exportfig(latmixroll, 'figures/latmixroll.eps', 'width', 32, 'color', 'cmyk','Fontmode','fixed','FontSize', 14);