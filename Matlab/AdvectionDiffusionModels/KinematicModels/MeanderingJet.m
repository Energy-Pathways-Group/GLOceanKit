classdef MeanderingJet < StreamfunctionModel
    %MeanderingJet Amy Bowers kinematic model of a meandering jet
    % We are going to use the "meandering jet" streamfunction that is in Amy
    % Bower's "A simple kinematic mechanism for mixing fluid parcels across a
    % meandering jet". See also Roger Samelson's "Fluid exchange across a
    % meandering jet".
    
    properties
        % Parameters from Bower 1991.
        L = 40e3; % width of the jet [m]
        U = 100*(1e3/86400); % maximum velocity of the jet [m/s]
        A = 50e3; % wave amplitude (amplitude in y-direction) [m]
        Lx = 350e3; % wave length of the jet (oscillation in x-direction) [m]
        cx = 0*10*(1e3/86400); % wave phase velocity of the jet [m/s]

% xlim = [0
% ylim
    end
    
    properties (Dependent)
        k
    end

    methods
        function self = MeanderingJet()
            self.xlim = [0 2*self.Lx];
            self.ylim = [-5*self.L 5*self.L];
            self.xIsPeriodic = 1;

            self.xVisLim = [0 2*self.Lx];
            self.yVisLim = [-5*self.L 5*self.L];
            
            self.name = 'Meandering jet';
        end

        function value = get.k(self)
            value = 2*pi/self.Lx;
        end

        function theta = theta(self,t,x)
            theta = self.k*(x-self.cx*t);
        end

        function gamma = gamma(self,t,x,y)
            gamma = (y-self.A*cos(self.theta(t,x))) / self.L ./ (1+(self.k*self.A*sin(self.theta(t,x))).^2).^(1/2);
        end

        function psi = psi(self,t,x,y)
            psi = self.U*self.L*(1-tanh(self.gamma(t,x,y)) );
        end

        function u = u(self,t,x,y)
            u = self.U * (1+(self.k*self.A*sin(self.theta(t,x))).^2).^(-1/2) .* sech(self.gamma(t,x,y)).^2;
        end
        
        function v = v(self,t,x,y)
            v = -self.U * (self.A*self.k*(1+self.A*self.A*self.k*self.k)*sin(self.theta(t,x))) .* (1+(self.k*self.A*sin(self.theta(t,x))).^2).^(-3/2) .* sech(self.gamma(t,x,y)).^2;
        end
    end
end

