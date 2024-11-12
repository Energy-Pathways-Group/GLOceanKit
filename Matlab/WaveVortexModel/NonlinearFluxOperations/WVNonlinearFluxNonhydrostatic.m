classdef WVNonlinearFluxNonhydrostatic < WVNonlinearFlux
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
        cos_alpha
        sin_alpha
        ApmD
        ApmW

        % AwD
        % AwZ
        % AwW
        % AwN
    end
    methods
        function self = WVNonlinearFluxNonhydrostatic(wvt,options)
            arguments
                wvt WVTransform {mustBeNonempty}
                options.uv_damp (1,1) double
                options.w_damp (1,1) double % characteristic speed used to set the damping. Try using wMax
                options.nu_xy (1,1) double
                options.nu_z (1,1) double
                options.r (1,1) double {mustBeNonnegative} = 0 % linear bottom friction, try 1/(200*86400) https://www.nemo-ocean.eu/doc/node70.html
                options.shouldUseBeta double {mustBeMember(options.shouldUseBeta,[0 1])} = 0
            end

            qgArgs = namedargs2cell(options);
            self@WVNonlinearFlux(wvt,qgArgs{:});

            k = shiftdim(wvt.k,-1);
            l = shiftdim(wvt.l,-1);
            kappa = sqrt(k.^2 + l.^2);
            self.cos_alpha = k./kappa;
            self.sin_alpha = l./kappa;
            self.cos_alpha(1) = 0;
            self.sin_alpha(1) = 0;

            signNorm = -2*(mod(wvt.j,2) == 1)+1; % equivalent to (-1)^j
            prefactor = signNorm * sqrt((wvt.g*wvt.Lz)/(2*(wvt.N0*wvt.N0 - wvt.f*wvt.f)));
            mj = (wvt.j*pi/wvt.Lz);
            self.ApmD = (mj/2) .* prefactor;
            self.ApmW = sqrt(-1) * (kappa/2) .* prefactor;

            % prefactor = signNorm * sqrt((wvt.g*wvt.Lz)/(2*(wvt.N0*wvt.N0 - wvt.f*wvt.f)));
            % self.AwD = prefactor .* (-sqrt(-1)*mj./(2*kappa));
            % self.AwZ = prefactor .* (-mj * wvt.f)./(2*kappa.*wvt.Omega);
            % self.AwW = prefactor .* (sqrt(-1)*kappa./2);
            % self.AwN = prefactor .* (-(wvt.N0*wvt.N0)*kappa./(2*wvt.Omega));
        end

        function [uNL,vNL,wNL,nNL] = spatialFlux(self,wvt)
            uNL = wvt.u .* wvt.diffX(wvt.u)   + wvt.v .* wvt.diffY(wvt.u)   + wvt.w .*  wvt.diffZF(wvt.u);
            vNL = wvt.u .* wvt.diffX(wvt.v)   + wvt.v .* wvt.diffY(wvt.v)   + wvt.w .*  wvt.diffZF(wvt.v);
            wNL = wvt.u .* wvt.diffX(wvt.w)   + wvt.v .* wvt.diffY(wvt.w)   + wvt.w .*  wvt.diffZG(wvt.w);
            nNL = wvt.u .* wvt.diffX(wvt.eta) + wvt.v .* wvt.diffY(wvt.eta) + wvt.w .* (wvt.diffZG(wvt.eta) + wvt.eta .* self.dLnN2);
        end

        function varargout = compute(self,wvt,varargin)

            [uNL,vNL,wNL,nNL] = self.spatialFlux(wvt);

            u_hat = wvt.transformFromSpatialDomainWithFourier(-uNL);
            v_hat = wvt.transformFromSpatialDomainWithFourier(-vNL);
            w_hat = wvt.transformFromSpatialDomainWithFourier(-wNL);
            n_hat = wvt.transformFromSpatialDomainWithFourier(-nNL);

            iK = sqrt(-1)*repmat(shiftdim(wvt.k,-1),wvt.Nz,1);
            iL = sqrt(-1)*repmat(shiftdim(wvt.l,-1),wvt.Nz,1);

            n_bar = wvt.transformFromSpatialDomainWithGg(n_hat);
            zeta_bar = wvt.transformFromSpatialDomainWithFg(iK .* v_hat - iL .* u_hat);
            A0NL = wvt.A0Z.*zeta_bar + wvt.A0N.*n_bar;

            delta_bar = self.ApmD .* (wvt.DCT * (self.cos_alpha .* u_hat + self.sin_alpha .* v_hat));
            w_bar = self.ApmW .* (wvt.DST * w_hat);
            nw_bar = wvt.transformWithG_wg(n_bar - A0NL);
            ApNL = delta_bar + w_bar + wvt.ApmN .* nw_bar;
            AmNL = delta_bar + w_bar - wvt.ApmN .* nw_bar;

            % Dbar = wvt.DCT*(iK .* u_hat + iL .* v_hat);
            % Zbar = wvt.DCT*(iK .* v_hat - iL .* u_hat);
            % Wbar = wvt.DST*(w_hat);
            % Nbar = wvt.DST*(n_hat);
            % Ap = self.AwD .* Dbar + self.AwZ .* Zbar + self.AwW .* Wbar + self.AwN .* Nbar;
            % Am = self.AwD .* Dbar - self.AwZ .* Zbar + self.AwW .* Wbar - self.AwN .* Nbar;

            ApNL(:,1) = wvt.transformFromSpatialDomainWithFio(u_hat(:,1) - sqrt(-1)*v_hat(:,1))/2;
            AmNL(:,1) = conj(ApNL(:,1));

            phase = exp(-wvt.iOmega*(wvt.t-wvt.t0));
            ApNL = self.damp .* wvt.Ap + ApNL .* phase;
            AmNL = self.damp .* wvt.Am + AmNL .* conj(phase);
            A0NL = ((self.damp + self.betaA0) .* wvt.A0 + A0NL);

            varargout = {ApNL,AmNL,A0NL};
        end
    end
end
