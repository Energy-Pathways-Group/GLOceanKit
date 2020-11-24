function [parameters,error] = EstimateLinearVelocityFieldParameters( x, y, t, parametersToEstimate)

shouldEstimateU0V0 = 0;
shouldEstimateUtVt = 0;
shouldEstimateStrain = 0;
shouldEstimateVorticity = 0;
shouldEstimateDivergence = 0;
nParameters = 0;

if isequal(class(parametersToEstimate),'ModelParameter')
    for i=1:length(parametersToEstimate)
       if  parametersToEstimate(i) == ModelParameter.u0v0
           shouldEstimateU0V0 = 1; nParameters = nParameters + 2;
       elseif  parametersToEstimate(i) == ModelParameter.utvt
           shouldEstimateUtVt = 1; nParameters = nParameters + 2;
       elseif  parametersToEstimate(i) == ModelParameter.strain
           shouldEstimateStrain = 1; nParameters = nParameters + 2;
       elseif  parametersToEstimate(i) == ModelParameter.vorticity
           shouldEstimateVorticity = 1; nParameters = nParameters + 1;
       elseif  parametersToEstimate(i) == ModelParameter.divergence
           shouldEstimateDivergence = 1; nParameters = nParameters + 1;
       end
    end
elseif isequal(class(parametersToEstimate),'char')
    if strcmp(parametersToEstimate,'strain-diffusive')
        shouldEstimateStrain = 1;
    elseif strcmp(parametersToEstimate,'vorticity-strain-diffusive')
        shouldEstimateStrain = 1;
        shouldEstimateVorticity = 1;
    end
end

u0 = 0;
v0 = 0;
ut = 0;
vt = 0;
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
tB = zeros(length(t)*nDrifters,1);
xB = zeros(length(t)*nDrifters,1);
yB = zeros(length(t)*nDrifters,1);
u = zeros(length(t)*nDrifters,1);
v = zeros(length(t)*nDrifters,1);
for iDrifter=1:nDrifters
    indices = (1:length(t)) + (iDrifter-1)*length(t);
    onesB(indices,:) = ones(size(x(:,iDrifter)));
    zerosB(indices,:) = zeros(size(x(:,iDrifter)));
    tB(indices,:) = t;
    xB(indices,:) = x(:,iDrifter);
    yB(indices,:) = y(:,iDrifter);
    u(indices) = dxdt(:,iDrifter);
    v(indices) = dydt(:,iDrifter);
end

Ru = []; Rv = [];
if shouldEstimateU0V0 == 1
    Ru = cat(2,Ru,onesB,zerosB);
    Rv = cat(2,Rv,zerosB,onesB);
end
if shouldEstimateUtVt == 1
    Ru = cat(2,Ru,tB,zerosB);
    Rv = cat(2,Rv,zerosB,tB);
end
if shouldEstimateStrain == 1
    Ru = cat(2,Ru,xB/2,yB/2);
    Rv = cat(2,Rv,-yB/2,xB/2);
end
if shouldEstimateVorticity == 1
    Ru = cat(2,Ru,-yB/2);
    Rv = cat(2,Rv,xB/2);
end
if shouldEstimateDivergence == 1
    Ru = cat(2,Ru,xB/2);
    Rv = cat(2,Rv,yB/2);
end

U = cat(1,u,v);
R = cat(1,Ru,Rv);
m = (R.' * R) \ (R.' * U);

p = 0;
if shouldEstimateU0V0 == 1
    p = p + 1; u0 = m(p);
    p = p + 1; v0 = m(p);
end
if shouldEstimateUtVt == 1
    p = p + 1; ut = m(p);
    p = p + 1; vt = m(p);
end
if shouldEstimateStrain == 1
    p = p + 1; sigma_n = m(p);
    p = p + 1; sigma_s = m(p);
end
if shouldEstimateVorticity == 1
    p = p + 1; zeta = m(p);
end
if shouldEstimateDivergence == 1
    p = p + 1; delta = m(p);
end

% 
% uv_model=R*m;
% u = uv_model(1:(length(uv_model)/2));
% v = uv_model((length(uv_model)/2+1):end);
% u = reshape(u,[length(t) nDrifters]);
% v = reshape(v,[length(t) nDrifters]);
% 
% u_res = dxdt - u;
% v_res = dydt - v;
% 
% u_bg = mean(u_res,2);
% v_bg = mean(v_res,2);
% 
% u_sm = u_res - u_bg;
% v_sm = v_res - v_bg;
% 
% x_sm = cumtrapz(t,u_sm);
% y_sm = cumtrapz(t,v_sm);



parameters.u0 = u0;
parameters.v0 = v0;
parameters.ut = ut;
parameters.vt = vt;
parameters.sigma_n = sigma_n;
parameters.sigma_s = sigma_s;
parameters.zeta = zeta;
parameters.delta = delta;
parameters.sigma = sqrt(sigma_n^2 + sigma_s^2);
parameters.theta = atan2(sigma_s,sigma_n)/2;
parameters.kappa = 0; %mean(x_sm(end,:).^2 + y_sm(end,:).^2)/(4*(t(end)-t(1)));
error = 0;

end

