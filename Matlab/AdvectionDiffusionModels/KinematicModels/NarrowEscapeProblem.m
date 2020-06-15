classdef NarrowEscapeProblem < KinematicModel
    %A box with an adjustable opening
    %   https://en.wikipedia.org/wiki/Narrow_escape_problem
    
    properties
        L = 10e3; % box length [m]
        delta = 1e3; % wall thickness---shouldn't matter
        W = 3e3; % escape window size
    end
    
    methods
        function self = NarrowEscapeProblem()
            
            self.xVisLim = 0.1*(self.L+2*self.delta)*[-1 1] + [-self.delta self.L+self.delta];
            self.yVisLim = 0.1*(self.L+2*self.delta)*[-1 1] + [-self.delta self.L+self.delta];
            self.visualScale = 1e3;
            
            self.RegenerateObstacles();
            
            self.name = 'Narrow escape';
        end
        
        function set.W(self,W)
            self.W = W;
            self.RegenerateObstacles();
        end
        
        function set.L(self,L)
            self.L = L;
            self.RegenerateObstacles();
        end
        
        function set.delta(self,delta)
            self.delta = delta;
            self.RegenerateObstacles();
        end
        
        function RegenerateObstacles(self)
            a = self.L;
            d = self.delta;
            wt = (a-self.W)/2 + self.W;
            wb = (a-self.W)/2;
            
            bottom = polyshape([-d -d a+d a+d], [0 -d -d 0]);
            right = polyshape([a a a+d a+d],[a+d -d -d a+d]);
            top = polyshape([-d -d a+d a+d],[a+d a a a+d]);
            lefttop = polyshape([-d -d 0 0],[a+d wt wt a+d]);
            leftbottom = polyshape([-d -d 0 0],[wb -d -d wb]);
            
            self.obstacles = union([bottom; right; top; lefttop; leftbottom]);
        end
        
    end
end