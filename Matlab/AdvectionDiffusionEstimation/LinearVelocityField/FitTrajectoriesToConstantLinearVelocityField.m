function [parameters,error] = FitTrajectoriesToConstantLinearVelocityField( x, y, t, model)

u0 = 0;
v0 = 0;
sigma_n = 0;
sigma_s = 0;
zeta = 0;
delta = 0;

nDrifters = size(x,2);

% Compute velocities with 2nd order accuracy
D = FiniteDifferenceMatrix(1,t,1,1,2);
dxdt = D*x;
dydt = D*y;

% Now put all the data together
onesB = zeros(length(t)*nDrifters,1);
zerosB = zeros(length(t)*nDrifters,1);
xB = zeros(length(t)*nDrifters,1);
yB = zeros(length(t)*nDrifters,1);
u = zeros(length(t)*nDrifters,1);
v = zeros(length(t)*nDrifters,1);
for iDrifter=1:nDrifters
    indices = (1:length(t)) + (iDrifter-1)*length(t);
    onesB(indices,:) = ones(size(x(:,iDrifter)));
    zerosB(indices,:) = zeros(size(x(:,iDrifter)));
    xB(indices,:) = x(:,iDrifter);
    yB(indices,:) = y(:,iDrifter);
    u(indices) = dxdt(:,iDrifter);
    v(indices) = dydt(:,iDrifter);
end

if strcmp(model,'strain-diffusive')
    Ru = cat(2,onesB,zerosB,xB/2,yB/2);
    Rv = cat(2,zerosB,onesB,-yB/2,xB/2);
    
    U = cat(1,u,v);
    R = cat(1,Ru,Rv);
    
    m = (R.' * R) \ (R.' * U);
    
    u0 = m(1);
    v0 = m(2);
    sigma_n = m(3);
    sigma_s = m(4);
    zeta = 0;
    delta = 0;
elseif strcmp(model,'vorticity-strain-diffusive')
    Ru = cat(2,onesB,zerosB,xB/2,yB/2,-yB/2);
    Rv = cat(2,zerosB,onesB,-yB/2,xB/2,xB/2);
    
    U = cat(1,u,v);
    R = cat(1,Ru,Rv);
    
    m = (R.' * R) \ (R.' * U);
    
    u0 = m(1);
    v0 = m(2);
    sigma_n = m(3);
    sigma_s = m(4);
    zeta = m(5);
    delta = 0;
end

uv_model=R*m;
u = uv_model(1:(length(uv_model)/2));
v = uv_model((length(uv_model)/2+1):end);
u = reshape(u,[length(t) nDrifters]);
v = reshape(v,[length(t) nDrifters]);

u_res = dxdt - u;
v_res = dydt - v;

u_bg = mean(u_res,2);
v_bg = mean(v_res,2);

u_sm = u_res - u_bg;
v_sm = v_res - v_bg;

x_sm = cumtrapz(t,u_sm);
y_sm = cumtrapz(t,v_sm);



parameters.u0 = u0;
parameters.v0 = v0;
parameters.sigma_n = sigma_n;
parameters.sigma_s = sigma_s;
parameters.zeta = zeta;
parameters.delta = delta;
parameters.sigma = sqrt(sigma_n^2 + sigma_s^2);
parameters.theta = atan2(sigma_s,sigma_n)/2;
parameters.kappa = mean(x_sm(end,:).^2 + y_sm(end,:).^2)/(4*(t(end)-t(1)));
error = 0;

end

