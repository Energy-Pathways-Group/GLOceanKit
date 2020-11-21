%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Jacknife
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma = 4e-6;
kappa = 0.5;
T = round(1/sigma/86400)*86400;
T = 3*86400;
dt = 3600;

velocityField = LinearVelocityField(sigma,0,0);
integrator = AdvectionDiffusionIntegrator(velocityField,kappa);

x = linspace(-500,500,5);
y = linspace(-500,500,5);
[x0,y0] = ndgrid(x,y);

x0 = [-500; -250; 0; 0; 0; 0; 0; 250; 500];
y0 = [0; 0; -500; -250; 0; 250; 500; 0; 0;];

x0 = x0+5000;
y0 = y0+5000;

totalIterations = 100;
kappaEst = zeros(totalIterations,2);
sigmaEst = zeros(totalIterations,2);
thetaEst = zeros(totalIterations,2);

for i=1:totalIterations
    [t,x,y] = integrator.particleTrajectories(x0,y0,T,dt);
        
    parameters = FitTrajectoriesToConstantLinearVelocityField( x-15000, y+1000, t, 'strain-diffusive' );
        
    kappaEst(i,1) = parameters.kappa;
    sigmaEst(i,1) = parameters.sigma;
    thetaEst(i,1) = parameters.theta;
    
    parameters = FitTrajectoriesToConstantLinearVelocityField( x, y, t, 'strain-diffusive' );
    
    kappaEst(i,2) = parameters.kappa;
    sigmaEst(i,2) = parameters.sigma;
    thetaEst(i,2) = parameters.theta;
end

tic
for iEst=1:2

    % Kernel density estimate of the distribution.
    % https://www.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation
    data = cat(2,sigmaEst(:,iEst).*cos(2*thetaEst(:,iEst)),sigmaEst(:,iEst).*sin(2*thetaEst(:,iEst)));
    [bandwidth,density,X,Y]=kde2d(data);
    
    pctTarget = flip(0.1:0.1:0.9);
    
    dLevels = DensityLevelForCDF(X,Y,density, pctTarget);
    
    pctLabels = cell(length(pctTarget),1);
    for i=1:length(pctTarget)
        pctLabels{i} = sprintf('%d %%',round(pctTarget(i)*100));
    end
    
    figure('Position',[50 50 1000 400])
    
    subplot(1,2,1)
    M = contourf(X,Y,density,dLevels);
    hold on
    for i=2:2:10
        rectangle('Position',[-i -i 2*i 2*i]*1e-6, 'Curvature', [1 1]);
    end
%     xlim([sigma_n_axis(1) sigma_n_axis(end)])
%     ylim([sigma_s_axis(1) sigma_s_axis(end)])
    
    xlim([-1e-5 1e-5])
    ylim([-1e-5 1e-5])
    
    axis equal
    xlabel('\sigma_n')
    ylabel('\sigma_s')
    title('Kernel density estimate')
    cb = colorbar('eastoutside');
    cb.Ticks = dLevels;
    cb.TickLabels = pctLabels;
    
    subplot(1,2,2)
    histogram((kappaEst(:,iEst)),10)
    xlabel('\kappa')
    title('diffusivity estimate')
    
end
toc

