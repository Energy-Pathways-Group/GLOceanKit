classdef TimeDependentLinearVelocityField < LinearVelocityField
    %LinearVelocityField Generate trajectories and solutions for a linear
    %velocity field
    
    properties
        t
    end
        
    methods
        function self = TimeDependentLinearVelocityField(t,zeta,sigma,theta,kappa,varargin)
            self.t = reshape(t,[],1);
            self.zeta = TimeDependentLinearVelocityField.NormalizeVariable(self.t,zeta);
            self.sigma = TimeDependentLinearVelocityField.NormalizeVariable(self.t,sigma);
            self.theta = TimeDependentLinearVelocityField.NormalizeVariable(self.t,theta);
            self.kappa = TimeDependentLinearVelocityField.NormalizeVariable(self.t,kappa);

            if length(varargin) == 2
                self.u0 = TimeDependentLinearVelocityField.NormalizeVariable(self.t,varargin{1});
                self.v0 = TimeDependentLinearVelocityField.NormalizeVariable(self.t,varargin{2});
            else
                error('Optional arguments (u0,v0) must both be given together');
            end
        end
    end

    methods (Static)
        function var = NormalizeVariable(t,var)
            if length(var) == 1
                var = var*ones(size(t));
            elseif length(var) == n
                var = reshape(var,[],1);
            else
                error('variable must be of length 1 or length(t).');
            end
        end
    end
end

