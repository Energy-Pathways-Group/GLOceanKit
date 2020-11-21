function [parameters,error] = FitTrajectoriesToConstantLinearVelocityField( x, y, t, model)


nDrifters = size(x,2);
x_com = 0*mean(x,2);
y_com = 0*mean(y,2);

% Rewrite the drifter positions in terms of the expansion point
q = x - x_com;
r = y - y_com;

% Compute velocities
dxdt = (x(2:end,:)-x(1:end-1,:))./diff(t);
dydt = (y(2:end,:)-y(1:end-1,:))./diff(t);
q = q(1:end-1,:) + diff(q)/2;
r = r(1:end-1,:) + diff(r)/2;

t = t(1:end-1);

% Now put all the data together
onesB = zeros(length(t)*nDrifters,1);
zerosB = zeros(length(t)*nDrifters,1);
xB = zeros(length(t)*nDrifters,1);
yB = zeros(length(t)*nDrifters,1);
u = zeros(length(t)*nDrifters,1);
v = zeros(length(t)*nDrifters,1);
for iDrifter=1:nDrifters
    indices = (1:length(t)) + (iDrifter-1)*length(t);
    onesB(indices,:) = ones(size(q(:,iDrifter)));
    zerosB(indices,:) = zeros(size(q(:,iDrifter)));
    xB(indices,:) = q(:,iDrifter);
    yB(indices,:) = r(:,iDrifter);
    u(indices) = dxdt(:,iDrifter);
    v(indices) = dydt(:,iDrifter);
end

if strcmp(model,'strain-diffusive')
    Ru = cat(2,onesB,zerosB,xB/2,yB/2);
    Rv = cat(2,zerosB,onesB,-yB/2,xB/2);
    
    U = cat(1,u,v);
    R = cat(1,Ru,Rv);
    
    m = (R.' * R) \ (R.' * U);
    
    ubg = m(1);
    vbg = m(2);
    sigma_n = m(3);
    sigma_s = m(4);
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



parameters.u0 = ubg;
parameters.v0 = vbg;
parameters.sigma = sqrt(sigma_n^2 + sigma_s^2);
parameters.theta = atan2(sigma_s,sigma_n)/2;
parameters.kappa = mean(x_sm(end,:).^2 + y_sm(end,:).^2)/(4*(t(end)-t(1)));
error = 0;

end

