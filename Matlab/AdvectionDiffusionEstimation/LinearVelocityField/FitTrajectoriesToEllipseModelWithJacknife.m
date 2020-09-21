function [parameters,error] = FitTrajectoriesToEllipseModelWithJacknife( x, y, t, model, k)

nDrifters = size(x,2);

if nargin < 5
    k = nDrifters-1;
end

permutations = nchoosek(1:nDrifters,k);
nPermutations = size(permutations,1);

kappa = zeros(nPermutations,1);
sigma = zeros(nPermutations,1);
zeta = zeros(nPermutations,1);

for iPermutation=1:nPermutations
    x_j = x(:,permutations(iPermutation,:));
    y_j = y(:,permutations(iPermutation,:));
    [~, ~, q, r] = CenterOfMass( x_j, y_j );
    [Mxx, Myy, Mxy] = SecondMomentMatrix( q, r);
    
    [parameters,error] = FitSecondMomentToEllipseModel( Mxx, Myy, Mxy, t, model);
    kappa(iPermutation) = parameters.kappa;
    sigma(iPermutation) = parameters.sigma*exp(2*sqrt(-1)*parameters.theta);
    zeta(iPermutation) = parameters.zeta;
end

parameters.kappa = mean(kappa);
parameters.sigma = abs(mean(sigma));
parameters.theta = angle(mean(sigma))/2;
parameters.zeta = mean(zeta);

end

