%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Unit testing the analytical path
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


nModels = 3;

for iModel = 1:nModels
    if iModel == 1
        modelName = 'strain';
        nParameters = 3;
        sigmaValues = [2e-6 8e-6 2e-5];
        zetaValues = 0*zeros(size(sigmaValues));
        thetaValues = [45 0 -45]*pi/180;
        T = 3*86400;
        dt = 3600;
    elseif iModel == 2
        modelName = 'vorticity-strain';
        nParameters = 3;
        sigmaValues = [2e-6 8e-6 2e-5];
        zetaValues = -0.5*sigmaValues;
        thetaValues = [45 0 -45]*pi/180;
        T = 5*86400;
        dt = 3600;
    elseif iModel == 3
        modelName = 'strain-vorticity';
        nParameters = 3;
        sigmaValues = [2e-6 8e-6 2e-5];
        zetaValues = -2*sigmaValues;
        thetaValues = [45 0 -45]*pi/180;
        T = 5*86400;
        dt = 3600;
    end
    
    for iParameter = 1:nParameters
        velocityField = LinearVelocityField(sigmaValues(iParameter),thetaValues(iParameter),zetaValues(iParameter));
        integrator = AdvectionDiffusionIntegrator(velocityField,0);
        
%         x = linspace(-500,500,15);
%         y = linspace(-100,100,15);
%         [x0,y0] = ndgrid(x,y);
        
         x0 = 4*[-500; -250; 0; 0; 0; 0; 0; 250; 500];
         y0 = [0; 0; -500; -250; 0; 250; 500; 0; 0;];
        
        [t,x,y] = integrator.particleTrajectories(x0,y0,T,dt/10);
        [x_a,y_a] = velocityField.ParticlePath(x0,y0,t,0,0,0);
        
        dx = x - x_a;
        dy = y - y_a;
        
        sqrt(mean(mean(dx.^2 + dy.^2)))
        
    end
end

