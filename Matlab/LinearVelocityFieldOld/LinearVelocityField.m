classdef LinearVelocityField < handle
    %LinearVelocityField Generate trajectories and solutions for a linear
    %velocity field
    
    properties
        sigma
        theta
        zeta
        kappa
        
        u0
        v0
    end
    
    properties (Dependent)
        sigma_n
        sigma_s
    end
    
    methods
        function self = LinearVelocityField(zeta,sigma,theta,kappa,varargin)
            self.zeta = zeta;
            self.sigma = sigma;
            self.theta = theta;
            self.kappa = kappa;
            
            if length(varargin) == 2
                self.u0 = varargin{1};
                self.v0 = varargin{2};
            elseif isempty(varargin) == 1
            else
                error('Optional arguments (u0,v0) must both be given together');
            end
        end
        
        function value = get.sigma_n(self)
            value = self.sigma*cos(2*self.theta);
        end
        
        function value = get.sigma_s(self)
            value = self.sigma*sin(2*self.theta);
        end
    end


    methods (Static)
        function [sigma_n,sigma_s] = NormalAndShearFromSigmaTheta(sigma,theta)
            sigma_n = sigma*cos(2*theta);
            sigma_s = sigma*sin(2*theta);
        end

        function [sigma,theta] = SigmaThetaFromNormalAndShear(sigma_n,sigma_s)
            sigma = sqrt(sigma_n.*sigma_n + sigma_s.*sigma_s);
            theta = atan2(sigma_n,sigma_s)/2;
        end
    end
end

