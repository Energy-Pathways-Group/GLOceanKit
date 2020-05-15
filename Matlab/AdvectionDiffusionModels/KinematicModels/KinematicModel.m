classdef KinematicModel < handle
    %KinematicModel Kinematic models are time-dependent 2-dimensional
    %velocity fields that do not require integration.
    
    properties
        xlim = [-Inf Inf] % box limits of the ocean/domain, Inf is valid
        ylim = [-Inf Inf] % box limits of the ocean/domain, Inf is valid
        xIsPeriodic = 0 % If xlim is finite, this may be set to 1 (true)
        yIsPeriodic = 0 % If ylim is finite, this may be set to 1 (true)
        obstacles = [] % polygon structure defining island/continents/other structures
        
        xVisLim % visual (box) limits of the ocean, no infinities
        yVisLim % visual (box) limits of the ocean, no infinities   
        
        name = ''
    end
    
    
    methods (Abstract)
       u = u(self,t,x,y);
       v = v(self,t,x,y);
    end

    methods
        function [x0,y0] = removeOutOfBoundsParticles(self,x0,y0)
            outOfBounds = x0 < min(self.xlim) | x0 > max(self.xlim) | y0 < min(self.ylim) | y0 > max(self.ylim);
            if ~isempty(self.obstacles)
                outOfBounds = outOfBounds | isinterior(self.obstacles,x0,y0);
            end
            if sum(outOfBounds) > 0
                fprintf('Removed %d particles because they were out of bounds.\n',sum(outOfBounds));
                x0(outOfBounds) = [];
                y0(outOfBounds) = [];
            end
        end
        
        function varargout = plotVelocityField(self,t)
            if ~exist('t','var')
                t = 0;
            end
            N = 15;
            xg = linspace(min(self.xVisLim),max(self.xVisLim),2*N)';
            yg = linspace(min(self.yVisLim),max(self.yVisLim),N)';
            [X,Y] = meshgrid(xg,yg);
                       
            if ~isempty(self.obstacles)
                x = reshape(X,[],1);
                y = reshape(Y,[],1);
                mask = ~isinterior(self.obstacles,x,y);
                mask = reshape(mask,size(X));

                quiver(X/1e3,Y/1e3,mask.*self.u(t,X,Y),mask.*self.v(t,X,Y))
                plot(scale(self.obstacles,1e-3))
            else
                quiver(X/1e3,Y/1e3,self.u(t,X,Y),self.v(t,X,Y))
            end
            
            axis equal
            xlim([min(self.xVisLim) max(self.xVisLim)]/1e3)
            ylim([min(self.yVisLim) max(self.yVisLim)]/1e3)
            xlabel('km'), ylabel('km')
            title(sprintf('%s t=%d',self.name,t))
            
            if length(nargout) == 2
                varargout{1} = X;
                varargout{2} = Y;
            end
        end
        
        function plotTrajectories(self,x,y)
            if self.xIsPeriodic == 1
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Here's a trick for drawing periodic data. If we use mod(x,xWindow), then
                % we will get weird wrapping effects
                xMin = min(self.xlim);
                xMax = max(self.xlim);
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

    methods (Static)
        function bounds = boundsFromLimits(xlim,ylim)
            xv = [min(xlim) min(xlim) max(xlim) max(xlim)];
            yv = [max(ylim) min(ylim) min(ylim) max(ylim)];
            bounds = struct('xv',xv,'yv',yv);
        end
    end
end

