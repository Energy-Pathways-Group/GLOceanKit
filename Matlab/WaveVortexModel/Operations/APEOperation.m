classdef APEOperation < WVOperation
% Computes the nonlinear flux for a WVTransform

    properties (GetAccess=public, SetAccess=protected)
        % rho_nm
        % p_nm
    end

    methods

        function self = APEOperation(wvt)
            arguments
                wvt WVTransform
            end
            outputVariables(1) = WVVariableAnnotation('ape',{'x','y','z'},'m^2 s^{-2}', 'available potential energy density');
            self@WVOperation('ape',outputVariables,@disp);

            % % In the definitions of rho_nm and p_nm we do not include the
            % % factor of rho0 that should be multiplying everything.
            % if isa(wvt.rhoFunction,'chebfun')
            %     self.rho_nm = wvt.rhoFunction/wvt.rho0 - 1;
            % else
            %     self.rho_nm = chebfun(wvt.rhoFunction,[min(wvt.z) max(wvt.z)],'splitting','on')/wvt.rho0 - 1;
            % end
            % self.p_nm = - wvt.g * cumsum(self.rho_nm); 
            % self.p_nm = self.p_nm - self.p_nm(0);
        end

        function varargout = compute(self,wvt,varargin)
            rho_nm = wvt.chebfunForZArray(wvt.rho_nm)/wvt.rho0 - 1;
            p_nm = - wvt.g * cumsum(rho_nm); 
            p_nm = p_nm - p_nm(0);

            eta_true = wvt.eta_true;            
            Z = wvt.Z;
            ape = wvt.g*eta_true.*rho_nm(Z - eta_true) + p_nm(Z) - p_nm(Z - eta_true);

            % This is an alternate way to compute ape, although it is much slower.
            % ape = zeros(wvt.spatialMatrixSize);
            % for i=1:numel(Z)
            %     z = Z(i);
            %     eta = eta_true(i);
            %     ape(i) = integral( @(xi) xi.*wvt.N2Function(z-xi),0,eta);
            % end

            varargout = {ape};
        end
  
    end

end