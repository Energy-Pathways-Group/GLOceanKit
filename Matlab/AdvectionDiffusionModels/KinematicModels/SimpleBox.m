classdef SimpleBox < KinematicModel
    %SimpleBox A box with no flow
    
    properties
        Lx = 10e3; % width of the jet [m]
        Ly = 5e3;
    end
    
    methods
        function self = SimpleBox()
            self.xlim = [0 self.Lx];
            self.ylim = [0 self.Ly];
            
            epsilonx = 0.05*self.Lx;
            epsilony = 0.05*self.Ly;
            self.xVisLim = [0-epsilonx self.Lx+epsilonx];
            self.yVisLim = [0-epsilony self.Ly+epsilony];
            
            self.name = 'No flow box';
        end
       
        function u = u(self,t,x,y)
            u = zeros(size(x));
        end
        
        function v = v(self,t,x,y)
            v = zeros(size(y));
        end
    end
end

