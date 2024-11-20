classdef WVHydrostaticFlux < WVNonlinearFluxOperation
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
        function self = WVHydrostaticFlux(wvt,force)
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            arguments (Repeating)
                force WVForcing
            end
            fluxVar(1) = WVVariableAnnotation('Fp',{'j','kl'},'m/s2', 'non-linear flux into Ap');
            fluxVar(2) = WVVariableAnnotation('Fm',{'j','kl'},'m/s2', 'non-linear flux into Am');
            fluxVar(3) = WVVariableAnnotation('F0',{'j','kl'},'m/s', 'non-linear flux into A0');

            self@WVNonlinearFluxOperation('WVHydrostaticFlux',fluxVar);

            if isa(wvt,'WVTransformConstantStratification')
                self.dLnN2 = 0;
            elseif isa(wvt,'WVTransformHydrostatic')
                self.dLnN2 = shiftdim(wvt.dLnN2,-2);
            elseif isa(wvt,'WVTransformBoussinesq')
                self.dLnN2 = shiftdim(wvt.dLnN2,-2);
            else
                self.dLnN2 = shiftdim(wvt.dLnN2,-2);
                warning('WVTransform not recognized.')
            end
            for iForce=1:length(force)
                self.addForcing(force{iForce});
            end
        end

        function varargout = compute(self,wvt,varargin)
            uNL = -(wvt.u .* wvt.diffX(wvt.u)   + wvt.v .* wvt.diffY(wvt.u)   + wvt.w .*  wvt.diffZF(wvt.u));
            vNL = -(wvt.u .* wvt.diffX(wvt.v)   + wvt.v .* wvt.diffY(wvt.v)   + wvt.w .*  wvt.diffZF(wvt.v));
            nNL = -(wvt.u .* wvt.diffX(wvt.eta) + wvt.v .* wvt.diffY(wvt.eta) + wvt.w .* (wvt.diffZG(wvt.eta) + wvt.eta .* self.dLnN2));

            for iForce=1:length(self.spatialForcing)
                [uNL,vNL,nNL] = self.spatialForcing{iForce}.addHydrostaticSpatialForcing(wvt, uNL, vNL, nNL);
            end
            
            [Fp,Fm,F0] = wvt.transformUVEtaToWaveVortex(uNL,vNL,nNL,wvt.t);

            for iForce=1:length(self.spectralForcing)
                [Fp,Fm,F0] = self.spectralForcing{iForce}.addSpectralForcing(wvt, Fp, Fm, F0);
            end
            varargout = {Fp,Fm,F0};
        end

        function addForcing(self,force)
            self.forcing{end+1} = force;
            if force.doesHydrostaticSpatialForcing == 1
                self.spatialForcing{end+1} = force;
            end
            if force.doesSpectralForcing == 1
                self.spectralForcing{end+1} = force;
            end
        end
    end

    methods (Static)
        function nlFlux = nonlinearFluxFromFile(ncfile,wvt)
            arguments
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            nlFlux = WVHydrostaticFlux(wvt);
            nlFlux.initForcingFromFile(ncfile,wvt);
        end
    end
end
