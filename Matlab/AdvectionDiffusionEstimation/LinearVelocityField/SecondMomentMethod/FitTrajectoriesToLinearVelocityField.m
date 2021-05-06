function [parameters,error] = FitTrajectoriesToLinearVelocityField( x, y, t, model, dof, K)

% First order of business is to set the expansion point.
t_knot = InterpolatingSpline.KnotPointsForSplines(t,K,dof);
nDrifters = size(x,2);
xf = zeros(size(x));
yf = zeros(size(x));
for iDrifter = 1:nDrifters
    xd = x(:,iDrifter);
    yd = y(:,iDrifter);

    spline_mean_x = ConstrainedSpline(t,xd,K,t_knot,NormalDistribution(1),[]);
    spline_mean_y = ConstrainedSpline(t,yd,K,t_knot,NormalDistribution(1),[]);
    
    xf(:,iDrifter) = spline_mean_x(t);
    yf(:,iDrifter) = spline_mean_y(t);
end
[x_com, y_com, ~, ~] = CenterOfMass(xf,yf);

% Rewrite the drifter positions in terms of the expansion point
q = x - x_com;
r = y - y_com;

% Compute velocities
dqdt = (q(2:end,:)-q(1:end-1,:))./diff(t);
drdt = (r(2:end,:)-r(1:end-1,:))./diff(t);
q = q(1:end-1) + diff(q)/2;
r = r(1:end-1) + diff(r)/2;

t = t(1:end-1);
t_knot = InterpolatingSpline.KnotPointsForSplines(t,K,dof);

% B is matrix of size(B)=[length(t) nSplines]
B = BSpline.Spline(t,t_knot,K,0);
nSplines = size(B,2);

% Now put all the data together
onesB = zeros(length(t)*nDrifters,nSplines);
zerosB = zeros(length(t)*nDrifters,nSplines);
xB = zeros(length(t)*nDrifters,nSplines);
yB = zeros(length(t)*nDrifters,nSplines);
u = zeros(length(t)*nDrifters,1);
v = zeros(length(t)*nDrifters,1);
for iDrifter=1:nDrifters
    indices = (1:length(t)) + (iDrifter-1)*length(t);
    onesB(indices,:) = ones(size(q(:,iDrifter))).*B;
    zerosB(indices,:) = zeros(size(q(:,iDrifter))).*B;
    xB(indices,:) = q(:,iDrifter).*B;
    yB(indices,:) = r(:,iDrifter).*B;
    u(indices) = dqdt;
    v(indices) = drdt;
end

if strcmp(model,'strain-diffusive')
    Ru = cat(2,onesB,zerosB,xB,yB);
    Rv = cat(2,zerosB,onesB,-yB,xB);
    
    U = cat(1,u,v);
    R = cat(1,Ru,Rv);
    
    m = (R.' * R) \ (R.' * U);
    
    ubg = B*m(1:nSplines);
    vbg = B*m( (1:nSplines)+nSplines);
    sigma_n = B*m( (1:nSplines)+nSplines);
    sigma_s = B*m( (1:nSplines)+nSplines);
end


end

