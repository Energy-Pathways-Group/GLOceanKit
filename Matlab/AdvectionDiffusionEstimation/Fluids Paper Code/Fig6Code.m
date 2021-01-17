clear all
load("smoothedGriddedRho1Drifters.mat") %using latmix initial positions as before
x = x(1,1:9); y = y(1,1:9);
NX = 9; % number of drifters
NT = 300; % number of time increments for simulation
dt = 1800; % time increment in seconds
kappa_sm=0.5;
%%
rng(456)
delta = 0;  
sigma = 1e-5:(1e-6-1e-5)/NT:(1e-6-(1e-6-1e-5)/NT); 
theta = 30;
zeta = 0; 
sigma_s = sigma*sind(2*theta);
sigma_n = sigma*cosd(2*theta);
for ii = 1:NT
    x(ii+1,:) = x(ii,:)+0.5*(sigma_n(ii)+delta)*x(ii,:)*dt+0.5*(sigma_s(ii)-zeta)*y(ii,:)*dt+sqrt(2*dt)*sqrt(kappa_sm)*randn(1,NX);
    y(ii+1,:) = y(ii,:)+0.5*(sigma_s(ii)+zeta)*x(ii,:)*dt+0.5*(delta-sigma_n(ii))*y(ii,:)*dt+sqrt(2*dt)*sqrt(kappa_sm)*randn(1,NX);
end
%%
u=diff(x)/dt; % calculate velocities
v=diff(y)/dt;
x = x(1:(size(x,1)-1),:); % we have to difference to get velocities, which loses a point - remove the final time point to make all variables the same size
y = y(1:(size(y,1)-1),:);
%%
NX2 = 100; NY2 = 100;
DX = max(max(abs(x))); DY = max(max(abs(y)));
[x2,y2]=meshgrid((DX/NX2)*(-(NX2-0.5)+(13*NX2)/12:NX2/12.5:NX2-0.5),(DY/NY2)*(-(NY2-0.5)+(13*NX2)/12:NX2/12.5:NY2-0.5));
x2 = x2/1000; y2=y2/1000;
u_mes = 0.5*(sigma_n(1)+delta).*x2+0.5*(sigma_s(1)-zeta).*y2;
v_mes = 0.5*(sigma_s(1)+zeta).*x2+0.5*(delta-sigma_n(1)).*y2;
fig1 = figure; set(gcf,'Position',[100 100 2100 600])
subplot(1,2,1);
quiver(x2,y2,u_mes,v_mes); hold on
plot(x/1000,y/1000)
xlabel('km'); ylabel('km'); axis equal
%%
x_com = x-mean(x,2); y_com = y-mean(y,2);
NX2 = 100; NY2 = 100;
DX = max(max(abs(x_com))); DY = max(max(abs(y_com)));
[x2,y2]=meshgrid((DX/NX2)*(-(NX2-0.5):NX2/12.5:NX2-0.5),(DY/NY2)*(-(NY2-0.5):NX2/12.5:NY2-0.5));
x2 = x2/1000; y2=y2/1000;
u_mes = 0.5*(sigma_n(1)+delta).*x2+0.5*(sigma_s(1)-zeta).*y2;
v_mes = 0.5*(sigma_s(1)+zeta).*x2+0.5*(delta-sigma_n(1)).*y2;
subplot(1,2,2); quiver(x2,y2,u_mes,v_mes); hold on
plot(x_com/1000, y_com/1000)
xlim([min(min(x_com/1000)) max(max(x_com/1000))]); ylim([min(min(y_com/1000)) max(max(y_com/1000))])
xlabel('km'); ylabel('km'); axis equal
exportfig(fig1, 'figures/PositionsVaryStrain.eps', 'width', 32, 'color', 'cmyk','Fontmode','fixed','FontSize', 14); 