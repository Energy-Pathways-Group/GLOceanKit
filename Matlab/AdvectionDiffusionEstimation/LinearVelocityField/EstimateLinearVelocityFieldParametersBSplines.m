function [parameters] = EstimateLinearVelocityFieldParametersBSplines( x, y, t, parametersToEstimate, B)
% EstimateLinearVelocityFieldParametersBSplines This function estimates
% linear velocity field parameters strain, vorticity and divergence from a
% cluster of Lagrangian particles and set of basis function B.
%
% Required inputs are,
%   x ? [nT nDrifters] x position in meters of the nDrifters
%   y - [nT nDrifters] y position in meters of the nDrifters
%   t - [nT 1] times in seconds of the observation times
%   parametersToEstimate - array of ModelParameter objects.
%   B - [nT nSplines] B-splines (or really any basis functions) to be used
%   for the parameters.
%
% Output is,
%   parameters - struct containing,
%       [u0,v0,u1,v1,sigma_n,sigma_s,zeta,delta] ? these will be [nT 1]
%       unless dof=1, in which case they will be scalar values.
%
% Note that this implements equations 18-20 in Oscroft, Sykulski & Early. A
% proper implementation would remove the artificial requirement that
% drifter observations occur at the same time.
%

shouldEstimateU0V0 = 0;
shouldEstimateU1V1 = 0;
shouldEstimateStrain = 0;
shouldEstimateVorticity = 0;
shouldEstimateDivergence = 0;
nParameters = 0;

K = max(sum(B > 0,2));

if isequal(class(parametersToEstimate),'ModelParameter')
    for i=1:length(parametersToEstimate)
       if  parametersToEstimate(i) == ModelParameter.u0v0
           shouldEstimateU0V0 = 1; nParameters = nParameters + 2;
       elseif  parametersToEstimate(i) == ModelParameter.u1v1
           if K > 1
              error('u1v1 is invalid with splines K>1!'); 
           end
           shouldEstimateU1V1 = 1; nParameters = nParameters + 2;
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

u0 = zeros(length(t),1);
v0 = zeros(length(t),1);
u1 = zeros(length(t),1);
v1 = zeros(length(t),1);
sigma_n = zeros(length(t),1);
sigma_s = zeros(length(t),1);
zeta = zeros(length(t),1);
delta = zeros(length(t),1);

nDrifters = size(x,2);

% Center-of-mass position, and positions relative to the center-of-mass
mx = mean(x,2);
my = mean(y,2);

q = x-mx;
r = y-my;
% Center-of-mass velocity, and velocities relative to the center-of-mass

if length(unique(diff(t))) > 1
    % Compute velocities with 2nd order accuracy
    D = FiniteDifferenceMatrix(1,t,1,1,2);
    
    dmxdt = D*mx;
    dmydt = D*my;
    
    dqdt = D*q;
    drdt = D*r;
else
    % 2nd order finite difference for evenly spaced dt
    dmxdt = zeros(length(t),1);
    dmydt = zeros(length(t),1);
    dqdt = zeros(length(t),nDrifters);
    drdt = zeros(length(t),nDrifters);
    dt = t(2)-t(1);
    
    dmiddle = @(x) (x(3:end,:)-x(1:end-2,:))/(2*dt);
    dleft = @(x) (-x(3,:) + 4*x(2,:) - 3*x(1,:))/(2*dt);
    dright = @(x) (x(end-2,:) - 4*x(end-1,:) + 3*x(end,:))/(2*dt);
    
    dmxdt(2:end-1) = dmiddle(mx);
    dmydt(2:end-1) = dmiddle(my);
    dqdt(2:end-1,:) = dmiddle(q);
    drdt(2:end-1,:) = dmiddle(r);
    
    dmxdt(1) = dleft(mx);
    dmydt(1) = dleft(my);
    dqdt(1,:) = dleft(q);
    drdt(1,:) = dleft(r);
    
    dmxdt(end) = dright(mx);
    dmydt(end) = dright(my);
    dqdt(end,:) = dright(q);
    drdt(end,:) = dright(r);
end

nSplines = size(B,2);

% Now put all the data together
onesB = zeros(length(t)*nDrifters,nSplines);
zerosB = zeros(length(t)*nDrifters,nSplines);
tB = zeros(length(t)*nDrifters,nSplines);
qB = zeros(length(t)*nDrifters,nSplines);
rB = zeros(length(t)*nDrifters,nSplines);
u = zeros(length(t)*nDrifters,1);
v = zeros(length(t)*nDrifters,1);
for iDrifter=1:nDrifters
    indices = (1:length(t)) + (iDrifter-1)*length(t);
    onesB(indices,:) = ones(size(q(:,iDrifter))).*B;
    zerosB(indices,:) = zeros(size(q(:,iDrifter))).*B;
    tB(indices,:) = t.*B;
    qB(indices,:) = q(:,iDrifter).*B;
    rB(indices,:) = r(:,iDrifter).*B;
    u(indices) = dqdt(:,iDrifter);
    v(indices) = drdt(:,iDrifter);
end

% needed for the COM trajectory.
onesM = ones(length(t),nSplines).*B;
zerosM = zeros(length(t),nSplines);
tM = t.*B;
mxB = mx.*B;
myB = my.*B;

% Ru/Rv stores the data for the drifters *in* the COM frame
% Ru_cm/Rv_cm stores the data *of* the COM frame.
Ru = []; Rv = []; Ru_cm = []; Rv_cm = [];
if shouldEstimateU0V0 == 1
    Ru = cat(2,Ru,zerosB,zerosB);
    Rv = cat(2,Rv,zerosB,zerosB);
    Ru_cm = cat(2,Ru_cm,onesM,zerosM);
    Rv_cm = cat(2,Rv_cm,zerosM,onesM);
end
if shouldEstimateU1V1 == 1
    Ru = cat(2,Ru,zerosB,zerosB);
    Rv = cat(2,Rv,zerosB,zerosB);
    Ru_cm = cat(2,Ru_cm,tM,zerosM);
    Rv_cm = cat(2,Rv_cm,zerosM,tM);
end
if shouldEstimateStrain == 1
    Ru = cat(2,Ru,qB/2,rB/2);
    Rv = cat(2,Rv,-rB/2,qB/2);
    Ru_cm = cat(2,Ru_cm,mxB/2,myB/2);
    Rv_cm = cat(2,Rv_cm,-myB/2,mxB/2);
end
if shouldEstimateVorticity == 1
    Ru = cat(2,Ru,-rB/2);
    Rv = cat(2,Rv,qB/2);
    Ru_cm = cat(2,Ru_cm,-myB/2);
    Rv_cm = cat(2,Rv_cm,mxB/2);
end
if shouldEstimateDivergence == 1
    Ru = cat(2,Ru,qB/2);
    Rv = cat(2,Rv,rB/2);
    Ru_cm = cat(2,Ru_cm,mxB/2);
    Rv_cm = cat(2,Rv_cm,myB/2);
end

U = cat(1,u,v);
R = cat(1,Ru,Rv);

U = cat(1,U,dmxdt,dmydt);
R = cat(1,R,Ru_cm,Rv_cm);

m = (R.' * R) \ (R.' * U);

p = 0;
if shouldEstimateU0V0 == 1
    p = p + 1; u0 = B*m((1:nSplines) + (p-1)*nSplines);
    p = p + 1; v0 = B*m((1:nSplines) + (p-1)*nSplines);
end
if shouldEstimateU1V1 == 1
    p = p + 1; u1 = B*m((1:nSplines) + (p-1)*nSplines);
    p = p + 1; v1 = B*m((1:nSplines) + (p-1)*nSplines);
end
if shouldEstimateStrain == 1
    p = p + 1; sigma_n = B*m((1:nSplines) + (p-1)*nSplines);
    p = p + 1; sigma_s = B*m((1:nSplines) + (p-1)*nSplines);
end
if shouldEstimateVorticity == 1
    p = p + 1; zeta = B*m((1:nSplines) + (p-1)*nSplines);
end
if shouldEstimateDivergence == 1
    p = p + 1; delta = B*m((1:nSplines) + (p-1)*nSplines);
end

if nSplines == 1
    % In this case there a a single spline, and it is just 1s. So return a
    % scalar value instead of the whole time series.
    parameters.u0 = u0(1);
    parameters.v0 = v0(1);
    parameters.u1 = u1(1);
    parameters.v1 = v1(1);
    parameters.sigma_n = sigma_n(1);
    parameters.sigma_s = sigma_s(1);
    parameters.zeta = zeta(1);
    parameters.delta = delta(1);
else
    parameters.u0 = u0;
    parameters.v0 = v0;
    parameters.u1 = u1;
    parameters.v1 = v1;
    parameters.sigma_n = sigma_n;
    parameters.sigma_s = sigma_s;
    parameters.zeta = zeta;
    parameters.delta = delta;
end
% parameters.sigma = sqrt(sigma_n.^2 + sigma_s.^2);
% parameters.theta = atan2(sigma_s,sigma_n)/2;
%parameters.kappa = 0; %mean(x_sm(end,:).^2 + y_sm(end,:).^2)/(4*(t(end)-t(1)));

end

