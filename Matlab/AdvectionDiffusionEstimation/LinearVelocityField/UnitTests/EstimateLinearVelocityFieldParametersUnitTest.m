%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 1. Generate particles in a linear velocity field
% 2. Fit estimate the linear velocity field parameters using 2nd moment
% 3. Plot the estimates for a range of ensembles
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nModels = 3;

for iModel = 1:nModels
    if iModel == 1
        % Three different parameters values in a strain only model.
        modelName = 'strain-diffusive';
        model = [ModelParameter.strain];
        nParameters = 3;
        kappaValues = [0.5 1 2];
        sigmaValues = [2e-6 8e-6 2e-5];
        zetaValues = 0*zeros(size(sigmaValues));
        thetaValues = [45 0 -45]*pi/180;
        T = 3*86400;
        dt = 4*3600;
    elseif iModel == 2
        % Parameters values for a velocity field with both strain and
        % vorticity, but strain dominated.
        modelName = 'vorticity-strain-diffusive';
        model = [ModelParameter.strain,ModelParameter.vorticity];
        nParameters = 3;
        kappaValues = [0.5 1 2];
        sigmaValues = [2e-6 8e-6 2e-5];
        zetaValues = -0.5*sigmaValues;
        thetaValues = 0*zeros(size(sigmaValues));
        T = 5*86400;
        dt = 4*3600;
    elseif iModel == 3
        % Parameters values for a velocity field with both strain and
        % vorticity, but vortity dominated.
        modelName = 'vorticity-strain-diffusive';
        model = [ModelParameter.strain,ModelParameter.vorticity];
        nParameters = 3;
        kappaValues = [0.5 1 2];
        sigmaValues = [2e-6 8e-6 2e-5];
        zetaValues = -2*sigmaValues;
        thetaValues = 0*zeros(size(sigmaValues));
        T = 5*86400;
        dt = 4*3600;
    end
    
    
    % u0 = 0.1*0;
    % v0 = -0.2*0;
    % sigma = 1e-5;
    % zeta = -3e-6;
    % theta = 30*pi/180;
    % kappa = 0.1;
    % T = round(1/sigma/86400)*86400;
    % T = 2*86400;
    % dt = 3600;
    % model='strain-diffusive';
    % model = [ModelParameter.strain];
    % % model='vorticity-strain-diffusive';
    
    for iParameter = 1:nParameters
        sigma = sigmaValues(iParameter);
        zeta = zetaValues(iParameter);
        kappa = kappaValues(iParameter);
        theta = thetaValues(iParameter);
        u0 = 0;
        v0 = 0;
        velocityField = LinearVelocityField(sigma,theta,zeta,u0,v0);
        integrator = AdvectionDiffusionIntegrator(velocityField,kappa);
        integrator.stepSize = 1800;
        
%         x = linspace(-500,500,5);
%         y = linspace(-500,500,5);
%         [x0,y0] = ndgrid(x,y);
        
        % Particles placed in a "plus" configuration.
        x0 = [-500; -250; 0; 0; 0; 0; 0; 250; 500];
        y0 = [0; 0; -500; -250; 0; 250; 500; 0; 0;];
        
        
        totalIterations = 100;
        u0Est = zeros(totalIterations,1);
        v0Est = zeros(totalIterations,1);
        sigma_nEst = zeros(totalIterations,1);
        sigma_sEst = zeros(totalIterations,1);
        kappaEst = zeros(totalIterations,1);
        zetaEst = zeros(totalIterations,1);
        
        for i=1:totalIterations
            [t,x,y] = integrator.particleTrajectories(x0,y0,T,dt);
            
            parameterEstimates = EstimateLinearVelocityFieldParameters( x, y, t, model );
            u0Est(i) = parameterEstimates.u0;
            v0Est(i) = parameterEstimates.v0;
            sigma_nEst(i) = parameterEstimates.sigma_n;
            sigma_sEst(i) = parameterEstimates.sigma_s;
            zetaEst(i) = parameterEstimates.zeta;
            
            [u_meso,v_meso,u_bg,v_bg,u_sm,v_sm] = DecomposeTrajectories(x, y, t, parameterEstimates);
            x_sm = cumtrapz(t,u_sm);
            y_sm = cumtrapz(t,v_sm);
            kappaEst(i) = mean(x_sm(end,:).^2 + y_sm(end,:).^2)/(4*(t(end)-t(1)));
        end
        
        
        % Kernel density estimate of the distribution.
        % https://www.mathworks.com/matlabcentral/fileexchange/17204-kernel-density-estimation
        data = cat(2,sigma_nEst,sigma_sEst);
        [bandwidth,density,X,Y]=kde2d(data);
        
        figure('Position', [50 50 900 400])
        subplot(2,4,[1 2 5 6])
        contourf(X,Y,density);
        hold on
        for i=2:2:10
            rectangle('Position',[-i -i 2*i 2*i]*1e-6, 'Curvature', [1 1]);
        end
        scatter(sigma_nEst, sigma_sEst, 5^2, 0*[1 1 1],'filled');
        scatter(sigma*cos(2*theta), sigma*sin(2*theta), 20^2, 1*[1 1 1],'LineWidth',5);
        scatter(sigma*cos(2*theta), sigma*sin(2*theta), 20^2, 0*[1 1 1],'LineWidth',3);
        axis equal
        xlabel('\sigma_n')
        ylabel('\sigma_s')
        
        % make the zeros be white
        cmap = colormap;
        cmap(1,:)=1;
        colormap(cmap)
        
        subplot(2,4,3)
        histogram(u0Est,'Normalization','pdf')
        vlines(u0,'k--')
        title('u_0')
        
        subplot(2,4,4)
        histogram(v0Est,'Normalization','pdf')
        vlines(v0,'k--')
        title('v_0')
        
        subplot(2,4,7)
        histogram(zetaEst,'Normalization','pdf')
        vlines(zeta,'k--')
        title('\zeta')
        
        subplot(2,4,8)
        histogram(kappaEst,'Normalization','pdf')
        vlines(kappa,'k--')
        title('\kappa')
    end
end
