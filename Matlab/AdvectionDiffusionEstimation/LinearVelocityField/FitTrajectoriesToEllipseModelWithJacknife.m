function [parameters,error] = FitTrajectoriesToEllipseModelWithJacknife( x, y, t, model)

nDrifters = size(x,2);

kappa = zeros(nDrifters,1);
sigma = zeros(nDrifters,1);
zeta = zeros(nDrifters,1);

for iDrifter=1:nDrifters
    x_j = x; x_j(:,iDrifter) = [];
    y_j = y; y_j(:,iDrifter) = [];
    [~, ~, q, r] = CenterOfMass( x_j, y_j );
    [Mxx, Myy, Mxy] = SecondMomentMatrix( q, r);
    
    [parameters,error] = FitSecondMomentToEllipseModel( Mxx, Myy, Mxy, t, model);
    if strcmp(model,'diffusive')
        kappa(iDrifter) = parameters.kappa;
    elseif strcmp(model,'strain-diffusive')
        kappa(iDrifter) = parameters.kappa;
        sigma(iDrifter) = parameters.sigma*exp(2*sqrt(-1)*parameters.theta);
    end
end

parameters.kappa = mean(kappa);
parameters.sigma = abs(mean(sigma));
parameters.theta = angle(mean(sigma))/2;

end

