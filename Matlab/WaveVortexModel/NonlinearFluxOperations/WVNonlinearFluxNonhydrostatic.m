classdef WVNonlinearFluxNonhydrostatic < WVNonlinearFluxOperation
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
        function self = WVNonlinearFluxNonhydrostatic(wvt)
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            fluxVar(1) = WVVariableAnnotation('Fp',{'j','kl'},'m/s2', 'non-linear flux into Ap');
            fluxVar(2) = WVVariableAnnotation('Fm',{'j','kl'},'m/s2', 'non-linear flux into Am');
            fluxVar(3) = WVVariableAnnotation('F0',{'j','kl'},'m/s', 'non-linear flux into A0');

            self@WVNonlinearFluxOperation('WVNonlinearFluxNonhydrostatic',fluxVar);

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
        end

        function varargout = compute(self,wvt,varargin)
            uNL = wvt.u .* wvt.diffX(wvt.u)   + wvt.v .* wvt.diffY(wvt.u)   + wvt.w .*  wvt.diffZF(wvt.u);
            vNL = wvt.u .* wvt.diffX(wvt.v)   + wvt.v .* wvt.diffY(wvt.v)   + wvt.w .*  wvt.diffZF(wvt.v);
            wNL = wvt.u .* wvt.diffX(wvt.w)   + wvt.v .* wvt.diffY(wvt.w)   + wvt.w .*  wvt.diffZG(wvt.w);
            nNL = wvt.u .* wvt.diffX(wvt.eta) + wvt.v .* wvt.diffY(wvt.eta) + wvt.w .* (wvt.diffZG(wvt.eta) + wvt.eta .* self.dLnN2);

            u_hat = wvt.transformFromSpatialDomainWithFourier(-uNL);
            v_hat = wvt.transformFromSpatialDomainWithFourier(-vNL);
            w_hat = wvt.transformFromSpatialDomainWithFourier(-wNL);
            n_hat = wvt.transformFromSpatialDomainWithFourier(-nNL);

            iK = sqrt(-1)*repmat(shiftdim(wvt.k,-1),wvt.Nz,1);
            iL = sqrt(-1)*repmat(shiftdim(wvt.l,-1),wvt.Nz,1);

            n_bar = wvt.transformFromSpatialDomainWithGg(n_hat);
            zeta_bar = wvt.transformFromSpatialDomainWithFg(iK .* v_hat - iL .* u_hat);
            A0 = wvt.A0Z.*zeta_bar + wvt.A0N.*n_bar;

            delta_bar = wvt.transformWithG_wg(wvt.h_0.*wvt.transformFromSpatialDomainWithFg(iK .* u_hat + iL .* v_hat));
            nw_bar = wvt.transformWithG_wg(n_bar - A0);
            Ap = wvt.ApmD .* delta_bar + wvt.ApmN .* nw_bar;
            Am = wvt.ApmD .* delta_bar - wvt.ApmN .* nw_bar;

            Ap(:,1) = wvt.transformFromSpatialDomainWithFio(u_hat(:,1) - sqrt(-1)*v_hat(:,1))/2;
            Am(:,1) = conj(Ap(:,1));

            phase = exp(-wvt.iOmega*(t-wvt.t0));
            Ap = Ap .* phase;
            Am = Am .* conj(phase);

            varargout = {Ap,Am,A0};
        end
    end
end