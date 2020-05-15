classdef StreamfunctionModel < KinematicModel
    %StreamfunctionModel Streamfunction models are time-dependent
    %2-dimensional velocity fields defined by a stream function
    
    methods (Abstract)
       psi = psi(self,t,x,y);
    end

    methods

    function varargout = plotStreamfunction(self,t)
        if ~exist('t','var')
            t = 0;
        end

        N = 150;
        xg = linspace(min(self.xVisLim),max(self.xVisLim),2*N)';
        yg = linspace(min(self.yVisLim),max(self.yVisLim),N)';
        [Xf,Yf] = meshgrid(xg,yg);

        if ~isempty(self.obstacles)
            x = reshape(Xf,[],1);
            y = reshape(Yf,[],1);
            maskf = ~isinterior(self.obstacles,x,y);
            maskf = reshape(maskf,size(Xf));
            a = Xf(maskf); b = Yf(maskf);
            psiValues = self.psi(t,a,b);
            levels = linspace(min(psiValues),max(psiValues),20);
            
            contour(Xf/1e3,Yf/1e3,maskf.*self.psi(t,Xf,Yf),levels), hold on
            plot(scale(self.obstacles,1e-3))
        else
            contour(Xf/1e3,Yf/1e3,self.psi(t,Xf,Yf))
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

    end
end

