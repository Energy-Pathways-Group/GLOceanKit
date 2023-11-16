classdef WVNonlinearFluxSpatial < WVNonlinearFluxOperation
    % 3D nonlinear flux for Boussinesq flow, computed in the spatial domain
    %
    % Computes the nonlinear flux for a Boussinesq model. This class is not
    % intended to be used for numerical modeling as it does not have any
    % antialiasing or damping, but is indended as an example. The
    % implementation is *simple* and follows directly from the equations of
    % motion, but it is not the fastest implementation. To compute
    % nonlinear fluxes appropriate for numerical modeling, use the
    % [WVNonlinearFlux](/classes/wvnonlinearflux/) class.
    %
    % - Topic: Initializing
    % - Declaration: WVNonlinearFluxSpatial < [WVNonlinearFluxOperation](/classes/wvnonlinearfluxoperation/)
    properties
        dLnN2 = 0
    end
    methods
        function self = WVNonlinearFluxSpatial(wvt)
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            fluxVar(1) = WVVariableAnnotation('Fp',{'k','l','j'},'m/s2', 'non-linear flux into Ap');
            fluxVar(2) = WVVariableAnnotation('Fm',{'k','l','j'},'m/s2', 'non-linear flux into Am');
            fluxVar(3) = WVVariableAnnotation('F0',{'k','l','j'},'m/s', 'non-linear flux into A0');

            self@WVNonlinearFluxOperation('WVNonlinearFluxSpatial',fluxVar);

            if isa(wvt,'WVTransformConstantStratification')
                self.dLnN2 = 0;
            elseif isa(wvt,'WVTransformHydrostatic')
                self.dLnN2 = wvt.dLnN2;
            else
                self.dLnN2 = shiftdim(wvt.dLnN2,-2);
                warning('WVTransform not recognized.')
            end
        end

        function varargout = compute(self,wvt,varargin)
            varargout = cell(1,self.nVarOut);
            uNL = wvt.u .* wvt.diffX(wvt.u)   + wvt.v .* wvt.diffY(wvt.u)   + wvt.w .*  wvt.diffZF(wvt.u);
            vNL = wvt.u .* wvt.diffX(wvt.v)   + wvt.v .* wvt.diffY(wvt.v)   + wvt.w .*  wvt.diffZF(wvt.v);
            nNL = wvt.u .* wvt.diffX(wvt.eta) + wvt.v .* wvt.diffY(wvt.eta) + wvt.w .* (wvt.diffZG(wvt.eta) + wvt.eta .* shiftdim(self.dLnN2,-2));

            [varargout{:}] = wvt.transformUVEtaToWaveVortex(-uNL,-vNL,-nNL,wvt.t);
        end
    end
end
