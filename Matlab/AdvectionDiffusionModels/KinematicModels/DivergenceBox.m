classdef DivergenceBox < KinematicModel
    %SimpleBox A box with no flow
    
    properties
        Lx = 10e3;
        Ly = 5e3;
        Lg = 0.5e3;
        n = 2
        m = 1
        U = 1.0
    end
    
    methods
        function self = DivergenceBox()
            self.xlim = [0 self.Lx];
            self.ylim = [0 self.Ly];
            
            epsilonx = 0.05*self.Lx;
            epsilony = 0.05*self.Ly;
            self.xVisLim = [0-epsilonx self.Lx+epsilonx];
            self.yVisLim = [0-epsilony self.Ly+epsilony];
            self.visualScale = 1e3;
            
            self.name = 'No flow box';
        end
       
        % phi = sin(n*pi*x/Lx)*sin(m*pi*y/Ly);
        
        
        function u = u(self,t,x,y)
            u = zeros(size(x));
            for i=1:self.n
                for j=1:self.m
                    xc = (2*i-1)*self.Lx/(self.n+2);
                    yc = j*self.Ly/(self.m+1);
                    r2 = (x-xc).^2 + (y-yc).^2;
                    gaussian= exp( -r2/(2*self.Lg*self.Lg) );
                    u = u + ((-1).^i) * exp(0.5)*self.U * ((x-xc)/self.Lg) .* gaussian;
                end
            end
        end
        
        function v = v(self,t,x,y)
            v = zeros(size(x));
            for i=1:self.n
                for j=1:self.m
                    xc = (2*i-1)*self.Lx/(self.n+2);
                    yc = j*self.Ly/(self.m+1);
                    r2 = (x-xc).^2 + (y-yc).^2;
                    gaussian= exp( -r2/(2*self.Lg*self.Lg) );
                    v = v + ((-1).^i) * exp(0.5)*self.U * ((y-yc)/self.Lg) .* gaussian;
                end
            end
        end
    end
end

