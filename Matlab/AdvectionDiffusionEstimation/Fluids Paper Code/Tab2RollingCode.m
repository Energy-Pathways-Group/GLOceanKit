clear all
load('smoothedGriddedRho1Drifters.mat')
ND = 9;
W = 48;
dt = t(2)-t(1);
x=x(:,1:ND)-mean(x(:,1:ND),2); y=y(:,1:ND)-mean(y(:,1:ND),2);
u = diff(x(:,1:ND))/dt; v = diff(y(:,1:ND))/dt;
x = x(1:length(x)-1,1:ND); y = y(1:length(y)-1,1:ND);
Models = [1 0 1 1; 0 1 1 1; 1 1 0 0; 1 0 0 0; 0 1 0 0; 0 0 0 0];
for ii = 1:size(Models,1)
    jj = W/2+1;
    [delta(ii,jj), zeta(ii,jj), sigman(ii,jj), sigmas(ii,jj), u1, v1, u_sm, v_sm] = DivVortStrainEst( x(jj-W/2:jj+W/2,:), y(jj-W/2:jj+W/2,:), u(jj-W/2:jj+W/2,:), v(jj-W/2:jj+W/2,:), dt, Models(ii,1), Models(ii,2), Models(ii,3), Models(ii,4));
    USM(1:jj,:) = u_sm(1:jj,:); VSM(1:jj,:) = v_sm(1:jj,:);
    for jj = W/2+2:length(t)-2-W/2
        [delta(ii,jj), zeta(ii,jj), sigman(ii,jj), sigmas(ii,jj), u1, v1, u_sm, v_sm] = DivVortStrainEst( x(jj-W/2:jj+W/2,:), y(jj-W/2:jj+W/2,:), u(jj-W/2:jj+W/2,:), v(jj-W/2:jj+W/2,:), dt, Models(ii,1), Models(ii,2), Models(ii,3), Models(ii,4));
        USM(jj,:) = u_sm(W/2+1,:); VSM(jj,:) = v_sm(W/2+1,:);
    end
    jj = length(t)-1-W/2;
    [delta(ii,jj), zeta(ii,jj), sigman(ii,jj), sigmas(ii,jj), u1, v1, u_sm, v_sm] = DivVortStrainEst( x(jj-W/2:jj+W/2,:), y(jj-W/2:jj+W/2,:), u(jj-W/2:jj+W/2,:), v(jj-W/2:jj+W/2,:), dt, Models(ii,1), Models(ii,2), Models(ii,3), Models(ii,4));
    USM(jj:length(t)-1,:) = u_sm(W/2+1:end,:); VSM(jj:length(t)-1,:) = v_sm(W/2+1:end,:);
    %%%
    FVU(ii) = mean(mean(USM.^2+VSM.^2,2))./mean(mean(u.^2+v.^2,2));
    diff1 = dt/(4*length(t)).*abs(sum(USM+1i*VSM,1)).^2;
    kappa(ii) = mean(diff1);
    FDU(ii) = mean(abs(sum(USM+1i*VSM,1)).^2)/mean(abs(sum(u+1i*v,1)).^2);
end

sigma = sqrt(sigman.^2+sigmas.^2);
theta = atan2d(sigmas,sigman)/2;
%%
 [kappa;  FVU; FDU]'
%%