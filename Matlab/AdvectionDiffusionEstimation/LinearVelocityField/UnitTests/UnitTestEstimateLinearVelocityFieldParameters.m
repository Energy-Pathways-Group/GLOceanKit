%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Unit testing the estimators
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nModels = 4;

for iModel = 2:2%1:nModels
    if iModel == 1
       modelName = 'diffusive';
       nParameters = 3;
       kappaValues = [0.5 1 2];
       sigmaValues = 0*zeros(size(kappaValues));
       zetaValues = 0*zeros(size(sigmaValues));
       thetaValues = 0*zeros(size(sigmaValues));
       T = 3*86400;
       dt = 3600;
    elseif iModel == 2
        modelName = 'strain-diffusive';
        nParameters = 3;
        kappaValues = [0.5 1 2];
        sigmaValues = [2e-6 8e-6 2e-5];
        zetaValues = 0*zeros(size(sigmaValues));
        thetaValues = [45 0 -45]*pi/180;
        T = 3*86400;
        dt = 4*3600;
    elseif iModel == 3
        modelName = 'vorticity-strain-diffusive';
        nParameters = 3;
        kappaValues = [0.5 1 2];
        sigmaValues = [2e-6 8e-6 2e-5];
        zetaValues = -0.5*sigmaValues;
        thetaValues = 0*zeros(size(sigmaValues));
        T = 5*86400;
        dt = 4*3600;
    elseif iModel == 4
        modelName = 'vorticity-strain-diffusive';
        nParameters = 3;
        kappaValues = [0.5 1 2];
        sigmaValues = [2e-6 8e-6 2e-5];
        zetaValues = -2*sigmaValues;
        thetaValues = 0*zeros(size(sigmaValues));
        T = 5*86400;
        dt = 4*3600;
    end
    
    for iParameter = 1:nParameters
        velocityField = LinearVelocityField(sigmaValues(iParameter),thetaValues(iParameter),zetaValues(iParameter));
        integrator = AdvectionDiffusionIntegrator(velocityField,kappaValues(iParameter));
        
%         x = linspace(-500,500,15);
%         y = linspace(-100,100,15);
%         [x0,y0] = ndgrid(x,y);
        
         x0 = 4*[-500; -250; 0; 0; 0; 0; 0; 250; 500];
         y0 = [0; 0; -500; -250; 0; 250; 500; 0; 0;];
        
        
        totalIterations = 100;
        kappaEst = zeros(totalIterations,1);
        sigmaEst = zeros(totalIterations,1);
        thetaEst = zeros(totalIterations,1);
        zetaEst = zeros(totalIterations,1);
        
        kappaAltEst = zeros(totalIterations,1);
        
        for i=1:totalIterations
            [t,x,y] = integrator.particleTrajectories(x0,y0,T,dt);
            
%             parameters = FitTrajectoriesToEllipseModel( x, y, t, modelName, 'searchAlgorithm', 'fminsearch' );
            parameters = FitTrajectoriesToLinearizedEllipseModel(x, y, t, modelName);
            
            kappaEst(i) = parameters.kappa;
            sigmaEst(i) = parameters.sigma;
            thetaEst(i) = mod(parameters.theta+pi/2,pi)-pi/2;
            zetaEst(i) = parameters.zeta;
            
            kappaAltEst(i) = ((mean(x(end,:).^2+y(end,:).^2))-(mean(x(1,:).^2+y(1,:).^2)))/(4*t(end));
        end
        
        fprintf('expected (kappa,sigma,theta,zeta) = (%.2f, %.2g, %.1f, %.2g), found (kappa,sigma,theta,zeta) = (%.2f±%.2f, %.2g±%.2g, %.1f±%.1f, %.2g±%.2g)\n',kappaValues(iParameter),sigmaValues(iParameter),thetaValues(iParameter), zetaValues(iParameter), mean(kappaEst), std(kappaEst), median(sigmaEst), std(sigmaEst), mean(thetaEst), std(thetaEst), median(zetaEst), std(zetaEst) );
%         fprintf('alt: %.2f\n',mean(kappaAltEst));
    end
end

return

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

