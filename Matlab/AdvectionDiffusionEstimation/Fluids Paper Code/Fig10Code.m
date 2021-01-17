clear all
load('smoothedGriddedRho1Drifters.mat')
%%
x=x(:,1:9);
y=y(:,1:9);
fig1 = figure;
subplot(2,2,1)
plot(x/1000,y/1000)
hold on;
 plot(mean(x(1,:))/1000,mean(y(1,:))/1000,'k*','MarkerSize',30)
 plot(mean(x(end,:))/1000,mean(y(end,:))/1000,'r*','MarkerSize',30)
xlabel('km'); ylabel('km')
title('Site 1')
axis equal
subplot(2,2,3)
plot(x/1000 - mean(x/1000,2), y/1000 - mean(y/1000,2))
xlabel('km'); ylabel('km')
axis equal
%%
load('smoothedGriddedRho2Drifters.mat')
subplot(2,2,2)
plot(x/1000,y/1000)
xlabel('km'); ylabel('km')
hold on;
 plot(mean(x(1,:))/1000,mean(y(1,:))/1000,'k*','MarkerSize',30)
 plot(mean(x(end,:))/1000,mean(y(end,:))/1000,'r*','MarkerSize',30)
title('Site 2')
axis equal
subplot(2,2,4)
plot(x/1000 - mean(x/1000,2), y/1000 - mean(y/1000,2))
xlabel('km'); ylabel('km')
axis equal
 ha=get(gcf,'children');
 set(ha(1),'position',[.5 .1 .32 .32])
 set(ha(2),'position',[.5 .5 .32 .32])
 set(ha(3),'position',[.1 .1 .32 .32])
 set(ha(4),'position',[.1 .5 .32 .32])
exportfig(fig1, 'figures/LatmixPositions.eps', 'width', 32, 'color', 'cmyk','Fontmode','fixed','FontSize', 14);