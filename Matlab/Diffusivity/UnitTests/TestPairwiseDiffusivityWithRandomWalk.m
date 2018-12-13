Nx = 10;
Ny = 10;
reps = Nx*Ny;
kappa = 1.0;

t = (0:1000)';
deltaT = t(2)-t(1);
N = length(t);
sigma = sqrt(2*kappa/deltaT);

dX = sigma*randn(N,reps);
dY = sigma*randn(N,reps);

rmsDisplacement = sqrt(kappa*max(t))

% If I don't put these in a grid, then I get a length scale dependent
% diffusivity. Why? Because I'm grabbing and sorting, biasing the sample.
dL = 10; % sets the artificial grid space

[x0, y0] = ndgrid( dL*(1:Nx), dL*(1:Ny) );
x0 = reshape(x0, 1, reps);
y0 = reshape(y0, 1, reps);

x = deltaT*cumsum(dX,1) + repmat(x0,N,1);
y = deltaT*cumsum(dY,1) + repmat(y0,N,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
figure
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[edges,Nedges] = CreateBinEdgesForInitialSeparation(t, x, y);

[r0,D2] = PairwiseMeanSquareSeparation( t, x, y, edges );
kappa_r = zeros(size(r0));
kappa_r_err = zeros(size(r0));
for iBin = 1:size(kappa_r,2)
    [slope, slope_err] = LinearLeastSquaresFit(t,D2(:,iBin),1);
    kappa_r(iBin) = slope/4;
    kappa_r_err(iBin) = slope_err/4;
end

kappa_r(isnan(r0)) = [];
kappa_r_err(isnan(r0)) = [];
r0(isnan(r0)) = [];

subplot(2,1,1)
plot(t/86400,D2/1e6)
xlabel('days')
ylabel('km^2')

title(sprintf('Diffusivity set to %.2g m^2/s. Relative diffusivity should be double.',kappa))

subplot(2,1,2)
errorbar(r0/1e3, kappa_r,2*kappa_r_err)
xlabel('km')
ylabel('m^2/s')
ylim([0 1.1*max(kappa_r+2*kappa_r_err)])

