classdef CylinderFlow < StreamfunctionModel
    %CylinderFlow Potential flow around a cylinder
    %   https://en.wikipedia.org/wiki/Potential_flow_around_a_circular_cylinder
    
    properties
        R = 60e3; % radius of the cylinder [m]
        U = 100*(1e3/86400); % maximum velocity of the cylinder [m/s]
    end
    
    methods
        function self = CylinderFlow()
            self.xlim = 3*self.R*[-1 1];
            self.ylim = 3*self.R*[-1 1];
            self.xIsPeriodic = 1;
            
            self.xVisLim = 3*self.R*[-1 1];
            self.yVisLim = 3*self.R*[-1 1];
            self.visualScale = 1e3;
            
            theta = linspace(0,2*pi-2*pi/30,30);
            self.obstacles = polyshape(self.R*cos(theta),self.R*sin(theta));
            
            self.name = 'Cylinder flow';
        end
        
        function psi = psi(self,t,x,y)
            % psi = U (r - R^2/r)*sin(theta)
            % psi = U (y - R^2*y/r^2)
            r2 = x.^2 + y.^2;
            psi = self.U * y .* ( 1 - self.R*self.R ./ r2 ) ;
        end

        function u = u(self,t,x,y)
            r2 = x.^2 + y.^2;
            u = self.U - self.U * self.R * self.R * (x.^2 - y.^2 ) ./ (r2.*r2);
        end
        
        function v = v(self,t,x,y)
            r2 = x.^2 + y.^2;
            v = - self.U * self.R * self.R * (2 * x .* y ) ./ (r2.*r2);
        end
    end
end