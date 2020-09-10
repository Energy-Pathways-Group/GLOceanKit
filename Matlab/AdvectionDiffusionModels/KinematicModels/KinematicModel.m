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
        visualScale = 1e3
        
        name = ''
        u_ = [];
        v_ = [];
    end
    
    
%     methods (Abstract)
%        u = u(self,t,x,y);
%        v = v(self,t,x,y);
%     end

    methods
        function u = u(self,t,x,y)
            if isempty(self.u_)
                u = zeros(size(x));
            else
                u = self.u_(t,x,y);
            end
        end
        
        function v = v(self,t,x,y)
            if isempty(self.v_)
                v = zeros(size(x));
            else
                v = self.v_(t,x,y);
            end
        end
        
        function [x0,y0] = removeOutOfBoundsParticles(self,x0,y0)
            outOfBounds = x0 < min(self.xlim) | x0 > max(self.xlim) | y0 < min(self.ylim) | y0 > max(self.ylim);
            if ~isempty(self.obstacles)
                if self.xIsPeriodic == 1
                    x_wrap = mod(x0-min(self.xlim),max(self.xlim)-min(self.xlim)) + min(self.xlim);
                else
                    x_wrap = x0;
                end
                if self.yIsPeriodic == 1
                    y_wrap = mod(y0-min(self.ylim),max(self.ylim)-min(self.ylim)) + min(self.ylim);
                else
                    y_wrap = y0;
                end
                for iObstacle=1:length(self.obstacles)
                    outOfBounds = outOfBounds | isinterior(self.obstacles(iObstacle),x_wrap,y_wrap);
                end
            end
            if sum(outOfBounds) > 0
                fprintf('Removed %d particles because they were out of bounds.\n',sum(outOfBounds));
                x0(outOfBounds) = [];
                y0(outOfBounds) = [];
            end
        end
        
        function plotBounds(self,varargin)
            if all(~isinf(self.xlim)) && all(~isinf(self.ylim)) && self.xIsPeriodic == 0 && self.yIsPeriodic == 0
                rectangle('Position',[min(self.xlim) min(self.ylim) max(self.xlim)-min(self.xlim) max(self.ylim)-min(self.ylim)]/self.visualScale, varargin{:});
            elseif all(~isinf(self.ylim)) && self.yIsPeriodic == 0
                x = [min(self.xVisLim) min(self.xVisLim); max(self.xVisLim) max(self.xVisLim)];
                y = [min(self.yVisLim) max(self.yVisLim); min(self.yVisLim) max(self.yVisLim)];
                line( x/self.visualScale, y/self.visualScale, varargin{:});
            end
        end
        
        function varargout = plotVelocityField(self,t,quiverscale,N)
            if ~exist('t','var')
                t = 0;
            end
            if ~exist('quiverscale','var')
                quiverscale = 1;
            end
            if ~exist('N','var')
                N = 150;
            end
            
            
            xg = linspace(min(self.xVisLim),max(self.xVisLim),2*N)';
            yg = linspace(min(self.yVisLim),max(self.yVisLim),N)';
            [X,Y] = meshgrid(xg,yg);
                       
            if ~isempty(self.obstacles)
%                 x = reshape(X,[],1);
%                 y = reshape(Y,[],1);
%                 mask = false(size(x));
%                 for iObstacle=1:length(self.obstacles)
%                     mask = mask & ~isinterior(self.obstacles(iObstacle),x,y);
%                 end
%                 mask = reshape(mask,size(X));

%                 quiver(X/self.visualScale,Y/self.visualScale,mask.*self.u(t,X,Y),mask.*self.v(t,X,Y)), hold on
                quiver(X/self.visualScale,Y/self.visualScale,self.u(t,X,Y),self.v(t,X,Y),quiverscale), hold on
                plot(scale(self.obstacles,1e-3))
            else
                quiver(X/self.visualScale,Y/self.visualScale,self.u(t,X,Y),self.v(t,X,Y),quiverscale)
            end
            
            axis equal
            xlim([min(self.xVisLim) max(self.xVisLim)]/self.visualScale)
            ylim([min(self.yVisLim) max(self.yVisLim)]/self.visualScale)
            xlabel('km'), ylabel('km')
%             title(sprintf('%s t=%d',self.name,t))
            
            if length(nargout) == 2
                varargout{1} = X;
                varargout{2} = Y;
            end
        end
        
        function plotTrajectories(self,x,y,varargin)
%             shouldShow
%             extraargs = {};
%             for k = 1:2:length(extraargs)
%                 if strcmp(extraargs{k}, 'method')
%                     self.method = extraargs{k+1};
%                     extraargs(k+1) = [];
%                     extraargs(k) = [];
%                     userSpecifiedMethod = 1;
%                     break;
%                 end
%             end
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
                    plot((xshift+xMin)/self.visualScale,yshift/self.visualScale,varargin{:})
                end
                
                nmin = floor(min(x(end,:)-xMin)/xWindowLength);
                nmax = floor(max(x(end,:)-xMin)/xWindowLength);
                for n=nmin:nmax
                    xshift = x(end,:) - xMin - n*xWindowLength;
                    yshift = y(end,:);
                    mask = xshift < 0 | xshift > xWindowLength;
                    xshift(mask) = nan;
                    yshift(mask) = nan;
                    
                    scatter( (xshift+xMin)/self.visualScale,yshift/self.visualScale, 8^2, 'k', 'fill')
                end
            else
                plot(x/self.visualScale,y/self.visualScale,varargin{:}), hold on
                scatter( x(end,:)/self.visualScale,y(end,:)/self.visualScale, 8^2, 'k', 'fill')
                axis equal
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

