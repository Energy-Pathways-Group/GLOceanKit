classdef LinearVelocityField < StreamfunctionModel
    %LinearVelocityField Generate trajectories and solutions for a linear
    %velocity field
    
    properties
        sigma
        theta
        zeta
        
        u0
        v0
    end
    
    properties (Dependent)
        sigma_n
        sigma_s
    end
    
    methods
        function self = LinearVelocityField(sigma,theta,zeta,varargin)
            self.zeta = zeta;
            self.sigma = sigma;
            self.theta = theta;
            
            self.xVisLim = 1e3*[-1 1];
            self.yVisLim = 1e3*[-1 1];
            
            if length(varargin) == 2
                self.u0 = varargin{1};
                self.v0 = varargin{2};
            elseif isempty(varargin) == 1
            else
                error('Optional arguments (u0,v0) must both be given together');
            end
        end
        
        function set.zeta(self,val)
            self.zeta = val;
            self.nameDidChange();
        end
        function set.sigma(self,val)
            self.sigma = val;
            self.nameDidChange();
        end
        function set.theta(self,val)
            self.theta = val;
            self.nameDidChange();
        end
        
        function value = nameDidChange(self)
            if self.sigma*self.sigma > self.zeta*self.zeta
                value = 'Strain dominanted linear velocity field';
            elseif self.sigma*self.sigma < self.zeta*self.zeta
                value = 'Vorticity dominanted linear velocity field';
            else
                value = 'Strain-vorticity matched linear velocity field';
            end
            self.name = value;
        end
        
        function value = get.sigma_n(self)
            value = self.sigma*cos(2*self.theta);
        end
        
        function value = get.sigma_s(self)
            value = self.sigma*sin(2*self.theta);
        end
        
        function psi = psi(self,t,x,y)
            % Pulled from Sarah Oscroft's thesis chapter
            psi = 0.25*(self.sigma_s + self.zeta)*x.*x - 0.5*self.sigma_n*x.*y - 0.25*(self.sigma_s - self.zeta)*y.*y;
        end
        
        function u = u(self,t,x,y)
            u = (1/2)*(self.sigma_n*x + (self.sigma_s - self.zeta)*y);
        end
        
        function v = v(self,t,x,y)
            v = (1/2)*((self.sigma_s + self.zeta)*x - self.sigma_n*y);
        end
        
        function [x, y] = ParticlePathInLinearVelocityField( x0, y0, t, kappa, u_0, v_0 )
            %% ParticlePathInLinearVelocityField
            %
            %	Jeffrey J. Early, 2014
            %
            % (x0, y0) is the initial position of the particle, in meters.
            % This may be a vector for multiple particles.
            % t is the time vector, in seconds.
            % zeta is the vorticity, in 1/seconds
            % sigma_n and sigma_s are the normal and shear strain rate, in 1/seconds
            % (u_0, v_0) are the mean velocity, in m/s
            % kappa is the diffusivity, in m^2/s
            
            % Force (x0, y0) to be row vectors [1 m]
            % Force t to be a column vector [n 1]
            x0 = reshape(x0,1,[]);
            y0 = reshape(y0,1,[]);
            t = reshape(t,[],1);
            
            % Determine the step size for the random walk necessary to produce kappa
            deltaT = t(2)-t(1);
            N = length(t);
            
            % Generate the random increments
            randAmp = sqrt(deltaT*2*kappa);
            dX = randAmp*randn(N,length(x0));
            dY = randAmp*randn(N,length(y0));
            
            sigma2 = self.sigma_n*self.sigma_n + self.sigma_s*self.sigma_s;
            s2 = sigma2 - self.zeta*self.zeta;
            s = sqrt(s2);
            
            x = ones( length(t), 1 )*x0;
            y = ones( length(t), 1 )*y0;
            
            for n=2:length(t)
                % row vector, [1 m]
                % the previous position
                xn = x(n-1,:);
                yn = y(n-1,:);
                
                if (sigma2>0 && s2 == 0)
                    % The homogenous part of the solution
                    x(n,:) = (1+self.sigma_n*deltaT/2)*xn + ((self.sigma_s-self.zeta)*deltaT/2)*yn;
                    y(n,:) = ((self.sigma_s+self.zeta)*deltaT/2)*xn + (1-self.sigma_n*deltaT/2)*yn;
                    
                elseif (s2 < 0)
                    s = sqrt(self.zeta*self.zeta-sigma2);
                    cos_t = cos(s*deltaT/2);
                    sin_t = sin(s*deltaT/2);
                    
                    % The homogenous part of the solution
                    x(n,:) = (cos_t + (self.sigma_n/s)*sin_t)*xn + ((self.sigma_s-self.zeta)/s)*sin_t*yn;
                    y(n,:) = ((self.sigma_s+self.zeta)/s)*sin_t*xn + (cos_t - (self.sigma_n/s)*sin_t)*yn;
                    
                    % The non-homogenous part of the solution
                    x(n,:) = x(n,:) + (2/s^2)*( (s*sin_t + self.sigma_n*(1-cos_t))*u_0 + (self.sigma_s - self.zeta)*(1-cos_t)*v_0);
                    y(n,:) = y(n,:) + (2/s^2)*( (self.sigma_s + self.zeta)*(1-cos_t)*u_0 + (s*sin_t - self.sigma_n*(1-cos_t))*v_0);
                    
                elseif (s2 > 0 )
                    % scalar, [1 1]
                    cosh_t = cosh(s*deltaT/2);
                    sinh_t = sinh(s*deltaT/2);
                    
                    % The homogenous part of the solution
                    x(n,:) = (cosh_t + (self.sigma_n/s)*sinh_t)*xn + ((self.sigma_s-self.zeta)/s)*sinh_t*yn;
                    y(n,:) = ((self.sigma_s+self.zeta)/s)*sinh_t*xn + (cosh_t - (self.sigma_n/s)*sinh_t)*yn;
                    
                    % The non-homogenous part of the solution
                    x(n,:) = x(n,:) + (2/s^2)*( (s*sinh_t + self.sigma_n*(cosh_t - 1))*u_0 + (self.sigma_s - self.zeta)*(cosh_t - 1)*v_0);
                    y(n,:) = y(n,:) + (2/s^2)*( (self.sigma_s + self.zeta)*(cosh_t - 1)*u_0 + (s*sinh_t - self.sigma_n*(cosh_t - 1))*v_0);
                else
                    x(n,:) = xn;
                    y(n,:) = yn;
                end
                
                % The diffusive part of the solution
                x(n,:) = x(n,:) + dX(n,:);
                y(n,:) = y(n,:) + dY(n,:);
            end
        end
        
        function [Mxx, Myy, Mxy] = MomentTensorEvolution(self, Mxx0, Myy0, Mxy0, t, kappa ) 
            % Force t to be a column vector [n 1]
            % if (isrow(t))
            % 	t=t';
            % end
            
            if self.sigma==0.0 && self.zeta == 0.0
                Mxx = 2*kappa*t + Mxx0;
                Myy = 2*kappa*t + Myy0;
                Mxy = 0*t + Mxy0;
            elseif self.zeta==0.0
                cos_t = cos(self.theta);
                sin_t = sin(self.theta);
                
                cos2 = cos_t*cos_t;
                sin2 = sin_t*sin_t;
                cossin = cos_t*sin_t;
                
                tks = 2*kappa/self.sigma;
                
                % First take the initial conditions and rotate them to the strain aligned coordinates.
                Mxx0_r = cos2*Mxx0 + sin2*Myy0 + 2*cos_t*sin_t*Mxy0;
                Myy0_r = sin2*Mxx0 + cos2*Myy0 - 2*cos_t*sin_t*Mxy0;
                Mxy0_r = -cossin*Mxx0 + cossin*Myy0 + (cos2 - sin2)*Mxy0;
                
                % Write down the solution in the rotated frame
                Maa = (Mxx0_r + tks)*exp(self.sigma*t) - tks;
                Mbb = (Myy0_r - tks)*exp(-self.sigma*t) + tks;
                Mab = Mxy0_r;
                
                % Finally rotate the solution back to the original reference frame.
                Mxx = cos2*Maa + sin2*Mbb - 2*cos_t*sin_t*Mab;
                Myy = sin2*Maa + cos2*Mbb + 2*cos_t*sin_t*Mab;
                Mxy = cossin*Maa - cossin*Mbb + (cos2 - sin2)*Mab;
            elseif self.sigma==0.0
                A = (Mxx0 + Myy0)/2;
                B = -Mxy0;
                C = (Mxx0 - Myy0)/2;
                
                Mxx = A + 2*kappa*t + B*sin(self.zeta*t) + C*cos(self.zeta*t);
                Myy = A + 2*kappa*t - B*sin(self.zeta*t) - C*cos(self.zeta*t);
                Mxy = -B*cos(self.zeta*t) + C*sin(self.zeta*t);
            elseif self.zeta*self.zeta < self.sigma*self.sigma
                cos_t = cos(self.theta);
                sin_t = sin(self.theta);
                
                cos2 = cos_t*cos_t;
                sin2 = sin_t*sin_t;
                cossin = cos_t*sin_t;
                
                s = sqrt( self.sigma*self.sigma - self.zeta*self.zeta);
                tks = 2*kappa*self.sigma/(s*s);
                
                % First take the initial conditions and rotate them to the strain aligned coordinates.
                Mxx1 = cos2*Mxx0 + sin2*Myy0 + 2*cos_t*sin_t*Mxy0;
                Myy1 = sin2*Mxx0 + cos2*Myy0 - 2*cos_t*sin_t*Mxy0;
                Mxy1 = -cossin*Mxx0 + cossin*Myy0 + (cos2 - sin2)*Mxy0;
                
                % Compute the coefficients to the solution in the rotate frame.
                A = (1+self.sigma/s)*Mxx1/2 - (self.zeta/s)*Mxy1 - (1-self.sigma/s)*Myy1/2 + tks;
                B = (1-self.sigma/s)*Mxx1/2 + (self.zeta/s)*Mxy1 - (1+self.sigma/s)*Myy1/2 + tks;
                C = -(self.zeta/s)*Mxx1 + (2*self.sigma/s)*Mxy1 - (self.zeta/s)*Myy1;
                
                % Now write down the solution in the rotated frame
                Maa = (A/2)*(1+self.sigma/s)*exp(s*t) + (B/2)*(1-self.sigma/s)*exp(-s*t) + (self.zeta/s)*C/2 - (2*kappa/(s*s))*(self.zeta*self.zeta*t + self.sigma);
                Mbb = -(A/2)*(1-self.sigma/s)*exp(s*t) - (B/2)*(1+self.sigma/s)*exp(-s*t) + (self.zeta/s)*C/2 - (2*kappa/(s*s))*(self.zeta*self.zeta*t - self.sigma);
                Mab = (A/2)*(self.zeta/s)*exp(s*t) - (B/2)*(self.zeta/s)*exp(-s*t) + (self.sigma/s)*C/2 - (2*kappa/(s*s))*(self.zeta*self.sigma*t);
                
                % And finally rotate the solution back to the original reference frame.
                Mxx = cos2*Maa + sin2*Mbb - 2*cos_t*sin_t*Mab;
                Myy = sin2*Maa + cos2*Mbb + 2*cos_t*sin_t*Mab;
                Mxy = cossin*Maa - cossin*Mbb + (cos2 - sin2)*Mab;
            elseif self.zeta*self.zeta == self.sigma*self.sigma
                cos_t = cos(self.theta);
                sin_t = sin(self.theta);
                
                cos2 = cos_t*cos_t;
                sin2 = sin_t*sin_t;
                cossin = cos_t*sin_t;
                
                % First take the initial conditions and rotate them to the strain aligned coordinates.
                Mxx1 = cos2*Mxx0 + sin2*Myy0 + 2*cos_t*sin_t*Mxy0;
                Myy1 = sin2*Mxx0 + cos2*Myy0 - 2*cos_t*sin_t*Mxy0;
                Mxy1 = -cossin*Mxx0 + cossin*Myy0 + (cos2 - sin2)*Mxy0;
                
                A = -2*self.zeta*Mxy1 + self.sigma*(Mxx1 + Myy1);
                B = Mxx1 - Myy1;
                C = Mxx1 + Myy1;
                
                % Now write down the solution in the rotated frame
                Maa = kappa*self.sigma^2*t.^3/3 + kappa*self.sigma*t.^2 + 2*kappa*t + A*(self.sigma*t.^2+2*t)/4 + B*(self.sigma*t+1)/2 + C/2;
                Mbb = kappa*self.sigma^2*t.^3/3 - kappa*self.sigma*t.^2 + 2*kappa*t + A*(self.sigma*t.^2-2*t)/4 + B*(self.sigma*t-1)/2 + C/2;
                Mab = kappa*self.sigma*self.zeta*t.^3/3 + A*(self.zeta*t.^2/4 - 1/(2*self.zeta)) + B*self.zeta*t/2 + self.sigma*C/(2*self.zeta);
                
                % And finally rotate the solution back to the original reference frame.
                Mxx = cos2*Maa + sin2*Mbb - 2*cos_t*sin_t*Mab;
                Myy = sin2*Maa + cos2*Mbb + 2*cos_t*sin_t*Mab;
                Mxy = cossin*Maa - cossin*Mbb + (cos2 - sin2)*Mab;
                
            elseif self.zeta*self.zeta > self.sigma*self.sigma
                cos_t = cos(self.theta);
                sin_t = sin(self.theta);
                
                cos2 = cos_t*cos_t;
                sin2 = sin_t*sin_t;
                cossin = cos_t*sin_t;
                
                s = sqrt( self.zeta*self.zeta - self.sigma*self.sigma);
                
                % First take the initial conditions and rotate them to the strain aligned coordinates.
                Mxx1 = cos2*Mxx0 + sin2*Myy0 + 2*cos_t*sin_t*Mxy0;
                Myy1 = sin2*Mxx0 + cos2*Myy0 - 2*cos_t*sin_t*Mxy0;
                Mxy1 = -cossin*Mxx0 + cossin*Myy0 + (cos2 - sin2)*Mxy0;
                
                A = Mxx1 - Myy1 - 4*kappa*self.sigma/s^2;
                B = (self.sigma/s)*Mxx1 - (2*self.zeta/s)*Mxy1 + (self.sigma/s)*Myy1;
                C = (self.zeta^2/s^2)*Mxx1 - (2*self.zeta*self.sigma/s^2)*Mxy1 + (self.zeta^2/s^2)*Myy1;
                
                Maa = (A/2)*(cos(s*t) + (self.sigma/s)*sin(s*t)) + (B/2)*(sin(s*t)-(self.sigma/s)*cos(s*t)) + C/2 + (2*kappa/s^2)*(self.sigma + self.zeta^2*t);
                Mbb = -(B/2)*(sin(s*t) + (self.sigma/s)*cos(s*t)) - (A/2)*(cos(s*t)-(self.sigma/s)*sin(s*t)) + C/2 + (2*kappa/s^2)*(-self.sigma + self.zeta^2*t);
                Mab = (A*self.zeta/(2*s))*sin(s*t) - (B*self.zeta/(2*s))*cos(s*t) + C*self.sigma/(2*self.zeta) + 2*kappa*self.sigma*self.zeta*t/s^2;
                
                % And finally rotate the solution back to the original reference frame.
                Mxx = cos2*Maa + sin2*Mbb - 2*cos_t*sin_t*Mab;
                Myy = sin2*Maa + cos2*Mbb + 2*cos_t*sin_t*Mab;
                Mxy = cossin*Maa - cossin*Mbb + (cos2 - sin2)*Mab;
            end
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

