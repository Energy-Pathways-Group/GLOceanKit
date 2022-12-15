classdef AdvectionDiffusionIntegrator
    %AdvectionDiffusionModel Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        kinematicModel
        kappa
        stepSize = 0;
    end
    
    methods
        function self = AdvectionDiffusionIntegrator(aKinematicModel,kappa)
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
            if self.kinematicModel.isSpherical == 0
                flux = @(t,p) cat(2,self.kinematicModel.u(t,p(:,1),p(:,2)),self.kinematicModel.v(t,p(:,1),p(:,2)));
            else
                flux = @(t,p) cat(2,self.kinematicModel.u(t,p(:,1),p(:,2))./cos(p(:,2)/self.kinematicModel.R),self.kinematicModel.v(t,p(:,1),p(:,2)));
            end

            [x0,y0] = self.kinematicModel.removeOutOfBoundsParticles(x0,y0);
            if isempty(x0)
                fprintf('There were no particles left to advect. Aborting.\n');
            end
            n = length(x0);
            
            p0 = cat(2,x0,y0);
            tn = round(T/dt);
            
            % For the integrator, a periodic domain should be treated as if
            % it's unbounded.
%             if self.kinematicModel.xIsPeriodic == 1
%                 xMin = -Inf;
%                 xMax = Inf;
%             else
                xMin = min(self.kinematicModel.xlim);
                xMax = max(self.kinematicModel.xlim);
%             end
%             if self.kinematicModel.yIsPeriodic == 1
%                 yMin = -Inf;
%                 yMax = Inf;
%             else
                yMin = min(self.kinematicModel.ylim);
                yMax = max(self.kinematicModel.ylim);
%             end
            
            if self.stepSize == 0
                singleIncrement = 1;
                step = dt;
            else
                singleIncrement = 0;
                step = self.stepSize;
            end

            integrator = IntegratorWithObstacles( flux, p0,step,self.kappa,[xMin yMin],[xMax yMax], self.kinematicModel.obstacles, [self.kinematicModel.xIsPeriodic self.kinematicModel.yIsPeriodic] );
            
            t = zeros(tn+1,1);
            x = zeros(length(t),n);
            y = zeros(length(t),n);
            
            x(1,:) = x0;
            y(1,:) = y0;
            for i=1:tn
                if singleIncrement == 1
                    integrator.StepForwardOneIncrement;
                else
                    integrator.StepForwardToTime(i*dt);
                end
                p = integrator.currentY;
                x(i+1,:)=p(:,1).';
                y(i+1,:)=p(:,2).';
                t(i+1)=integrator.currentTime;
            end
        end

        function [t,x,y] = particleTrajectoriesWithStochasticUV(self,x0,y0,T,dt,u,v)
            % u and v are functions that takes (t) and returns a vector [n
            % 1] containing the velocity (u,v) for each particle. n is the
            % length of x0 and y0.

            % Returned values have the following shape,
            % size(t) = [nT 1], size(x) = [nT nParticles]
            x0 = reshape(x0,[],1);
            y0 = reshape(y0,[],1);
            if self.kinematicModel.isSpherical == 0
                flux = @(t,p) cat(2,u(t) + self.kinematicModel.u(t,p(:,1),p(:,2)),v(t) + self.kinematicModel.v(t,p(:,1),p(:,2)));
            else
                flux = @(t,p)  cat(2,(u(t) + self.kinematicModel.u(t,p(:,1),p(:,2)))./cos(p(:,2)/self.kinematicModel.R),v(t) + self.kinematicModel.v(t,p(:,1),p(:,2)));
            end

            [x0,y0] = self.kinematicModel.removeOutOfBoundsParticles(x0,y0);
            if isempty(x0)
                fprintf('There were no particles left to advect. Aborting.\n');
            end
            n = length(x0);
            
            p0 = cat(2,x0,y0);
            tn = round(T/dt);
            
            % For the integrator, a periodic domain should be treated as if
            % it's unbounded.
            xMin = min(self.kinematicModel.xlim);
            xMax = max(self.kinematicModel.xlim);
            yMin = min(self.kinematicModel.ylim);
            yMax = max(self.kinematicModel.ylim);

            if self.stepSize == 0
                singleIncrement = 1;
                step = dt;
            else
                singleIncrement = 0;
                step = self.stepSize;
            end

            integrator = IntegratorWithObstacles( flux, p0,step,self.kappa,[xMin yMin],[xMax yMax], self.kinematicModel.obstacles, [self.kinematicModel.xIsPeriodic self.kinematicModel.yIsPeriodic] );
            
            t = zeros(tn+1,1);
            x = zeros(length(t),n);
            y = zeros(length(t),n);
            
            x(1,:) = x0;
            y(1,:) = y0;
            for i=1:tn
                if singleIncrement == 1
                    integrator.StepForwardOneIncrement;
                else
                    integrator.StepForwardToTime(i*dt);
                end
                p = integrator.currentY;
                x(i+1,:)=p(:,1).';
                y(i+1,:)=p(:,2).';
                t(i+1)=integrator.currentTime;
            end
        end

    end
end

