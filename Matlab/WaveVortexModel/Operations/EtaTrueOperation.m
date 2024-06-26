classdef EtaTrueOperation < WVOperation
% Computes the nonlinear flux for a WVTransform

    properties (GetAccess=public, SetAccess=protected)
    
    end

    methods

        function self = EtaTrueOperation()
            outputVariables(1) = WVVariableAnnotation('eta_true',{'x','y','z'},'m', 'true isopycnal deviation');
            self@WVOperation('eta_true',outputVariables,@disp);
        end

        function varargout = compute(self,wvt,varargin)
            rho_nm = @(z) -(wvt.rhoFunction(z) - wvt.rho0);
            rho_total = -rho_nm(wvt.Z) + wvt.rho_e;


            zMinusEta = EtaTrueOperation.fInverseBisection(rho_nm,-rho_total(:),-wvt.Lz,0,1e-12);
            zMinusEta = reshape(zMinusEta,size(wvt.X));
            eta_true = wvt.Z - zMinusEta;
            varargout = {eta_true};
        end
  
    end

    methods (Static)
        function y = fInverseBisection(f, x, yMin,yMax, tol)
            %FINVERSEBISECTION(F, X)   Compute F^{-1}(X) using Bisection.
            % Taken from cumsum as part of chebfun.
            % chebfun/inv.m
            %
            % Copyright 2017 by The University of Oxford and The Chebfun Developers.
            % See http://www.chebfun.org/ for Chebfun information.

            a = yMin*ones(size(x));
            b = yMax*ones(size(x));
            c = (a + b)/2;

            while ( norm(b - a, inf) >= tol )
                vals = feval(f, c);
                % Bisection:
                I1 = ((vals-x) <= -tol);
                I2 = ((vals-x) >= tol);
                I3 = ~I1 & ~I2;
                a = I1.*c + I2.*a + I3.*c;
                b = I1.*b + I2.*c + I3.*c;
                c = (a+b)/2;
            end

            y = c;

        end
    end

end