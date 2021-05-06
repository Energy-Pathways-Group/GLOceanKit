function [parameters,error] = FitSecondMomentToEllipseModelLeastSquares( Mxx, Myy, Mxy, t, model)

% K = 3; % order of spline
% dSpline = 2*86400;
% 
% t_days = cat(1,(0:dSpline:t(end)).',t(end)); % observation points
% t_days(end-1) = []; % get rid of this point
% 
% t_knot = InterpolatingSpline.KnotPointsForPoints(t_days,K);
% 
% % B is matrix of size(B)=[length(t) nSplines]
% B = BSpline.Spline(t,t_knot,K,0);
% nSplines = size(B,2);

Mxx_t = vdiff(t(2)-t(1),Mxx,1);
Myy_t = vdiff(t(2)-t(1),Myy,1);
Mxy_t = vdiff(t(2)-t(1),Mxy,1);
nT = length(t);

kappaScale = 10^round(log10(mean( (Mxx_t+Myy_t) )));
lengthScale = 10^round(log10(sqrt(Mxx(end)+Myy(end))));
timeScale = lengthScale^2/kappaScale;
Mxx_t = Mxx_t/kappaScale;
Myy_t = Myy_t/kappaScale;
Mxy_t = Mxy_t/kappaScale;
Mxx = Mxx/(lengthScale^2);
Myy = Myy/(lengthScale^2);
Mxy = Mxy/(lengthScale^2);

% dt = t(2)-t(1);
% t_interp = t(1:(end-1)) + dt/2;
% Mxx_t = diff(Mxx)/dt;
% Myy_t = diff(Myy)/dt;
% Mxy_t = diff(Mxy)/dt;
% 
% nT = length(t_interp);

% By setting minError to Inf, these *will* be replaced.
if strcmp(model,'diffusive')
    % for [mxx + myy, mxx-myy, mxy]
    Xxx = 4*ones(nT,1);
    Xyy = zeros(nT,1);
    Xxy = zeros(nT,1);
    
    % for [mxx, myy, mxy]
    Xxx = 2*ones(nT,1);
    Xyy = 2*ones(nT,1);
    Xxy = zeros(nT,1);
elseif strcmp(model,'strain-diffusive')
    % for [mxx + myy, mxx-myy, mxy]
    Xxx = cat(2, 4*ones(nT,1),  Mxx-Myy,      2*Mxy);
    Xyy = cat(2, zeros(nT,1),   Mxx+Myy,      zeros(nT,1));
    Xxy = cat(2, zeros(nT,1),   zeros(nT,1),  0.5*(Mxx+Myy));
    
    % for [mxx, myy, mxy]
    Xxx = cat(2, 2*ones(nT,1),  Mxx,          Mxy);
    Xyy = cat(2, 2*ones(nT,1), -Myy,          Mxy);
    Xxy = cat(2, zeros(nT,1),   zeros(nT,1),  0.5*(Mxx+Myy));
    
elseif strcmp(model,'vorticity-strain*-diffusive') % strain dominanted => strain rate must be greater than zeta

elseif strcmp(model,'vorticity*-strain-diffusive') % strain dominanted => strain rate must be greater than zeta
    
end

R = cat(1,Xxx,Xyy,Xxy);
U = cat(1,Mxx_t,Myy_t,Mxy_t);
% U = cat(1,Mxx_t+Myy_t,Mxx_t-Myy_t,Mxy_t);
m = (R.' * R) \ (R.' * U);

if strcmp(model,'diffusive')
    parameters.kappa = m(1)*kappaScale;
    parameters.sigma = 0;
    parameters.theta = 0;
elseif strcmp(model,'strain-diffusive')
    parameters.kappa = m(1)*kappaScale;
    parameters.sigma = sqrt( m(2)^2 + m(3)^2 )/timeScale;
    parameters.theta = atan2(m(3),m(2))/2;
elseif strcmp(model,'vorticity-strain*-diffusive') % strain dominanted => strain rate must be greater than zeta
    
elseif strcmp(model,'vorticity*-strain-diffusive') % strain dominanted => strain rate must be greater than zeta
    
end

% D2 = (Mxx + Myy)/2;
% kappa_endpoint = 0.5 *( D2(end) - D2(1) )/ ( t(end) - t(1) )

error = mean((U - R*m).^2);

end

% kappa_raw = mean(diff(Mxx+Myy))/dt;
% t_diff = t(1:end-1)+dt/2;
% t_sum = cumsum(ones(size(diff(Mxx))))*dt;
% figure
% plot(t,Mxx+Myy), hold on
% plot(t_diff,cumsum(diff(Mxx+Myy)/dt)*dt)
% plot(t_diff,kappa_raw*t_sum)