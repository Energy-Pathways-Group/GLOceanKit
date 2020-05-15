classdef TranslatingGaussian < StreamfunctionModel
    %TranslatingGaussian A steady-translating Gaussian streamfunction.
    %   Simple eddy model.
    
    properties
        L = 60e3; % width of the jet [m]
        U = 0.12; % maximum velocity of the eddy [m/s]
        cx = -.0267; % translation velocity of the eddy [m/s]
        cy = 0; % translation velocity of the eddy [m/s]
    end
    
    methods
        function self = TranslatingGaussian()
            self.xVisLim = 3*self.L*[-1 1];
            self.yVisLim = 3*self.L*[-1 1];
            
            self.name = 'Translating Gaussian';
        end
        
        function psi = psi(self,t,x,y)
            % The factor of exp(0.5) fixes the amplitude so that max velocity is U
            r2 = (x-self.cx*t).^2 + (y-self.cy*t).^2;
            gaussian= exp( -r2/(2*self.L*self.L) );
            psi = exp(0.5)*self.U*self.L*gaussian;
        end

        function u = u(self,t,x,y)
            r2 = (x-self.cx*t).^2 + (y-self.cy*t).^2;
            gaussian= exp( -r2/(2*self.L*self.L) );
            u = exp(0.5)*self.U * ((y-self.cy*t)/self.L) .* gaussian;
        end
        
        function v = v(self,t,x,y)
            r2 = (x-self.cx*t).^2 + (y-self.cy*t).^2;
            gaussian= exp( -r2/(2*self.L*self.L) );
            v = -exp(0.5)*self.U * ((x-self.cx*t)/self.L) .* gaussian;
        end
    end
end

