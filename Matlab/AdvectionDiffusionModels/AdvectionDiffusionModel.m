classdef AdvectionDiffusionModel
    %AdvectionDiffusionModel Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        kinematicModel
        kappa
    end
    
    methods
        function self = AdvectionDiffusionModel(aKinematicModel,kappa)
            self.kinematicModel = aKinematicModel;
            self.kappa = kappa;
        end
        
        function [t,x,y] = particleTrajectories(self,x0,y0,T,dt)
            % (x0,y0) initial positions. Out-of-bounds particles will be
            % eliminated.
            % T is the total time the simulation is to be run
            % dt is the time step
            % 
            % Returned values have the following shape,
            % size(t) = [nT 1], size(x) = [nT nParticles]
            x0 = reshape(x0,[],1);
            y0 = reshape(y0,[],1);
            flux = @(t,p) cat(2,self.kinematicModel.u(t,p(:,1),p(:,2)),self.kinematicModel.v(t,p(:,1),p(:,2)));

            [x0,y0] = self.kinematicModel.removeOutOfBoundsParticles(x0,y0);
            if isempty(x0)
                fprintf('There were no particles left to advect. Aborting.\n');
            end
            n = length(x0);
            
            p0 = cat(2,x0,y0);
            tn = round(T/dt);
            
            % For the integrator, a periodic domain should be treated as if
            % it's unbounded.
            if self.kinematicModel.xIsPeriodic == 1
                xMin = -Inf;
                xMax = Inf;
            else
                xMin = min(self.kinematicModel.xlim);
                xMax = max(self.kinematicModel.xlim);
            end
            if self.kinematicModel.yIsPeriodic == 1
                yMin = -Inf;
                yMax = Inf;
            else
                yMin = min(self.kinematicModel.ylim);
                yMax = max(self.kinematicModel.ylim);
            end
            
            integrator = IntegratorWithDiffusivity( flux, p0,dt,self.kappa,[xMin yMin],[xMax yMax] );
            
            t = zeros(tn,1);
            x = zeros(length(t),n);
            y = zeros(length(t),n);
            
            x(1,:) = x0;
            y(1,:) = y0;
            for i=1:tn
                integrator.StepForwardToTime(i*integrator.stepSize);
                p = integrator.currentY;
                x(i,:)=p(:,1).';
                y(i,:)=p(:,2).';
                t(i)=integrator.currentTime;
            end
        end
        
        function plotTrajectories(self,x,y)
            if self.kinematicModel.xIsPeriodic == 1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Here's a trick for drawing periodic data. If we use mod(x,xWindow), then
                % we will get weird wrapping effects
                xMin = min(self.kinematicModel.xlim);
                xMax = max(self.kinematicModel.xlim);
                xWindowLength = xMax - xMin;
                nmin = floor(min(x(:)-xMin)/xWindowLength);
                nmax = floor(max(x(:)-xMin)/xWindowLength);
                h = gca;
                for n=nmin:nmax
                    xshift = x - xMin - n*xWindowLength;
                    yshift = y;
                    mask = xshift < 0 | xshift > xWindowLength;
                    xshift(mask) = nan;
                    yshift(mask) = nan;
                    h.ColorOrderIndex = 1;
                    plot((xshift+xMin)/1e3,yshift/1e3)
                end
            else
                plot(x/1e3,y/1e3)
            end
        end

    end
end

