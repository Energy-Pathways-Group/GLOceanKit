%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1. Generate particles in a linear velocity field
% 2. Fit estimate the linear velocity field parameters using 2nd moment
% 3. Plot the estimates for a range of ensembles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sigma = 2e-6;
theta = 0*pi/180;
kappa = 1.0;
T = round(1/sigma/86400)*86400;
T = 2*86400;
dt = 3600;

velocityField = LinearVelocityField(sigma,theta,0);
integrator = AdvectionDiffusionIntegrator(velocityField,kappa);

x = linspace(-500,500,5);
y = linspace(-500,500,5);
[x0,y0] = ndgrid(x,y);

x0 = [-500; -250; 0; 0; 0; 0; 0; 250; 500];
y0 = [0; 0; -500; -250; 0; 250; 500; 0; 0;];


totalIterations = 100;
kappaEst = zeros(totalIterations,1);
sigmaEst = zeros(totalIterations,1);
thetaEst = zeros(totalIterations,1);

for i=1:totalIterations
    [t,x,y] = integrator.particleTrajectories(x0,y0,T,dt);
        
    parameters = FitTrajectoriesToEllipseModel( x, y, t, 'strain-diffusive' );
        
    kappaEst(i) = parameters.kappa;
    sigmaEst(i) = parameters.sigma;
    thetaEst(i) = parameters.theta;
end

sigma_n = sigmaEst.*cos(2*thetaEst);
sigma_s = sigmaEst.*sin(2*thetaEst);

% Kernel density estimate of the distribution.
% https://www.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation
data = cat(2,sigma_n,sigma_s);
[bandwidth,density,X,Y]=kde2d(data);

figure
contourf(X,Y,density);
hold on
for i=2:2:10
    rectangle('Position',[-i -i 2*i 2*i]*1e-6, 'Curvature', [1 1]);
end
scatter(sigma_n, sigma_s, 5^2, 0*[1 1 1],'filled');
axis equal
xlabel('\sigma_n')
ylabel('\sigma_s')

% make the zeros be white
cmap = colormap;
cmap(1,:)=1;
colormap(cmap)

    return
    
    % Now let's identify the cumulative values enclosed by the different
    % contours.
    nLevels = 25;
    level = zeros(nLevels,1);
    pctEnclosed = zeros(nLevels,1);
    sigma_n_axis = X(1,:).';
    sigma_s_axis = Y(:,1);
    M = contourc(sigma_n_axis,sigma_s_axis,density,nLevels);
    i = 1;
    iLevel = 1;
    while (i < size(M,2))
        level(iLevel) = M(1,i);
        n = M(2,i);
        
        % the last point is redudant, and polyshape doesn't like that.
        cont = polyshape(M(1,(i+1):(i+n-1)),M(2,(i+1):(i+n-1)));
        
        mask = isinterior(cont,reshape(X,[],1),reshape(Y,[],1));
        mask = reshape(mask,size(X));
        
        pctEnclosed(iLevel) = trapz(sigma_s_axis,trapz(sigma_n_axis,density.*mask));
        
        i = i+n+1;
        iLevel = iLevel + 1;
    end
    
    pctTarget = flip(0.1:0.1:0.9);
    dLevels = interp1(pctEnclosed,level,pctTarget);
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
    xlim([sigma_n_axis(1) sigma_n_axis(end)])
    ylim([sigma_s_axis(1) sigma_s_axis(end)])
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

