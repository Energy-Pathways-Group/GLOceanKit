classdef EtaTrueOperation < WVOperation
% Computes the true vertical displacement, eta.
%
% 2025-07-14
% Note: this will recover rho_e with
%   rho_e = wvt.rhoFunction(wvt.Z-wvt.eta_true)-wvt.rhoFunction(wvt.Z);
% for which
%   drho = wvt.rho_e - rho_e;
%   sqrt(mean(drho(:).^2))
% gives 1e-13 for a realistic simulation.
%
% Using
%   int_vol = @(integrand) sum(mean(mean(shiftdim(wvt.z_int,-2).*integrand,1),2),3);
% I find that
%   int_vol(eta_true)
% returns ~4000 (or an average of 1 m displacement).
%
% This should be compared with
%   int_vol(shiftdim(wvt.N2,-2).*wvt.eta)/wvt.N2(end)
% which returns ~-65 or O(100) m.
%
% This all assume an adiabatic simulation---this is really just a measure
% of the MDA.

    properties (GetAccess=public, SetAccess=protected)
    
    end

    methods

        function self = EtaTrueOperation()
            outputVariables(1) = WVVariableAnnotation('eta_true',{'x','y','z'},'m', 'true isopycnal deviation');
            self@WVOperation('eta_true',outputVariables,@disp);
        end

        function varargout = compute(self,wvt,varargin)
            
            % rho_total = -rho_nm(wvt.Z) + wvt.rho_e;

            % Alternative 1: reference to the original no-motion state
            % rho_nm = @(z) -(wvt.rhoFunction(z) - wvt.rho0);

            % Alternative 2: reference to the current no-motion state
            % rho_nm_v = wvt.rho_nm - wvt.rho_nm(end);
            % rho_bar = squeeze(mean(mean(wvt.rho_e,1),2));
            % rho_nm_t = sort(rho_bar+rho_nm_v,1,"descend") - rho_nm_v;
            rho_nm_t = wvt.rho_nm - wvt.rho_nm0;
            rho_nm = @(z) -(wvt.rhoFunction(z) - wvt.rho0) - interp1(wvt.z,rho_nm_t,z,"linear");

            rho_total = (wvt.rhoFunction(wvt.Z) - wvt.rho0) + wvt.rho_e ;

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