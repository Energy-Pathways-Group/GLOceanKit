clear all
load('smoothedGriddedRho1Drifters.mat')
ND =9;
f0 = 2*(7.2921e-5)*sin(lat0*pi/180);
dt = t(2)-t(1);
x=x(:,1:ND)-mean(x(:,1:ND),2); y=y(:,1:ND)-mean(y(:,1:ND),2);
u = diff(x(:,1:ND))/dt; v = diff(y(:,1:ND))/dt;
x = x(1:length(x)-1,1:ND); y = y(1:length(y)-1,1:ND);
Models = [1 0 1 1; 0 1 1 1; 1 1 0 0; 1 0 0 0; 0 1 0 0; 0 0 0 0];
for ii = 1:size(Models,1)
 [delta(ii), zeta(ii), sigman(ii), sigmas(ii), u1, v1, u_sm, v_sm] = DivVortStrainEst( x, y, u, v, dt, Models(ii,1), Models(ii,2), Models(ii,3), Models(ii,4));
 FVU(ii) = mean(mean(u_sm.^2+v_sm.^2,2))./mean(mean(u.^2+v.^2,2));
 diff1 = dt/(4*length(t)).*abs(sum(u_sm+1i*v_sm,1)).^2;
 kappa(ii) = mean(diff1);
 FDU(ii) = mean(abs(sum(u_sm+1i*v_sm,1)).^2)/mean(abs(sum(u+1i*v,1)).^2);
end
sigma = sqrt(sigman.^2+sigmas.^2);
theta = atan2d(sigmas,sigman)/2;
[sigma/f0; theta; zeta/f0; delta/f0; kappa; FVU; FDU]'