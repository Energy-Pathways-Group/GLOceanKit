classdef Boussinesq < WVNonlinearFluxOperation

    properties
        shouldAntialias = 0
        AA
        nu_xy
        nu_z
        damp
    end

    methods
        function self = Boussinesq(wvt,options)
            arguments
                wvt WVTransform {mustBeNonempty}
                options.uv_damp (1,1) double % characteristic speed used to set the damping. Try using uMax.
                options.w_damp (1,1) double % characteristic speed used to set the damping. Try using wMax
                options.nu_xy (1,1) double
                options.nu_z (1,1) double
                options.shouldAntialias double = 0
                options.shouldUseSpectralVanishingViscosity = 1
            end
            fluxVar(1) = WVVariableAnnotation('Fp',{'k','l','j'},'m/s2', 'non-linear flux into Ap');
            fluxVar(2) = WVVariableAnnotation('Fm',{'k','l','j'},'m/s2', 'non-linear flux into Am');
            fluxVar(3) = WVVariableAnnotation('F0',{'k','l','j'},'m/s', 'non-linear flux into A0');

            self@WVNonlinearFluxOperation('Boussinesq',fluxVar);
            self.shouldAntialias = options.shouldAntialias;
            
            if self.shouldAntialias == 1
                self.AA = ~(wvt.maskForAliasedModes(jFraction=2/3));
                wvt.Ap = self.AA .* wvt.Ap;
                wvt.Am = self.AA .* wvt.Am;
                wvt.A0 = self.AA .* wvt.A0;
            else
                self.AA = 1;
            end
            

            if isfield(options,'nu_xy')
                self.nu_xy = nu_xy;
            else
                if isfield(options,'uv_damp')
                    if self.shouldAntialias == 1
                        self.nu_xy = (3/2)*(wvt.x(2)-wvt.x(1))*options.uv_damp/(pi^2);
                    else
                        self.nu_xy = (wvt.x(2)-wvt.x(1))*options.uv_damp/(pi^2);
                    end
                else
                    self.nu_xy = 0;
                end
            end

            if isfield(options,'nu_z')
                self.nu_z = nu_z;
            else
                if isfield(options,'w_damp')
                    if self.shouldAntialias == 1
                        self.nu_z = (3/2)*(wvt.z(2)-wvt.z(1))*options.w_damp/(pi^2);
                    else
                        self.nu_z = (wvt.z(2)-wvt.z(1))*options.w_damp/(pi^2);
                    end
                else
                    self.nu_z = 0;
                end
            end
            
            [K,L,J] = ndgrid(wvt.k,wvt.l,wvt.j);
            M = J*pi/wvt.Lz;
            self.damp = -(self.nu_z*M.^2 + self.nu_xy*(K.^2 +L.^2));
            if options.shouldUseSpectralVanishingViscosity == 1
                Qkl = wvt.spectralVanishingViscosityFilter(shouldAssumeAntialiasing=self.shouldAntialias);
                self.damp = Qkl.*self.damp;
            end

        end

        function varargout = compute(self,wvt,varargin)
            phase = exp(wvt.iOmega*(wvt.t-wvt.t0));
            Apt = wvt.Ap .* phase;
            Amt = wvt.Am .* conj(phase);
            A0t = wvt.A0;

            % Apply operator S---defined in (C4) in the manuscript
            Ubar = wvt.UAp.*Apt + wvt.UAm.*Amt + wvt.UA0.*A0t;
            Vbar = wvt.VAp.*Apt + wvt.VAm.*Amt + wvt.VA0.*A0t;
            Wbar = wvt.WAp.*Apt + wvt.WAm.*Amt;
            Nbar = wvt.NAp.*Apt + wvt.NAm.*Amt + wvt.NA0.*A0t;

            % Finishing applying S, but also compute derivatives at the
            % same time
            [U,Ux,Uy,Uz] = wvt.transformToSpatialDomainWithFAllDerivatives(Ubar);
            [V,Vx,Vy,Vz] = wvt.transformToSpatialDomainWithFAllDerivatives(Vbar);
            W = wvt.transformToSpatialDomainWithG(Wbar);
            [ETA,ETAx,ETAy,ETAz] = wvt.transformToSpatialDomainWithGAllDerivatives(Nbar);

            % compute the nonlinear terms in the spatial domain
            % (pseudospectral!)
            uNL = -U.*Ux - V.*Uy - W.*Uz;
            vNL = -U.*Vx - V.*Vy - W.*Vz;
            nNL = -U.*ETAx - V.*ETAy - W.*(ETAz + ETA.*shiftdim(wvt.dLnN2,-2));

            % Now apply the operator S^{-1} and then T_\omega^{-1}
            uNLbar = wvt.transformFromSpatialDomainWithF(uNL);
            vNLbar = wvt.transformFromSpatialDomainWithF(vNL);
            nNLbar = wvt.transformFromSpatialDomainWithG(nNL);

            Fp = self.AA .* (self.damp .* wvt.Ap + (wvt.ApU.*uNLbar + wvt.ApV.*vNLbar + wvt.ApN.*nNLbar) .* conj(phase));
            Fm = self.AA .* (self.damp .* wvt.Am + (wvt.AmU.*uNLbar + wvt.AmV.*vNLbar + wvt.AmN.*nNLbar) .* phase);
            F0 = self.AA .* (self.damp .* wvt.A0 + (wvt.A0U.*uNLbar + wvt.A0V.*vNLbar + wvt.A0N.*nNLbar));

            varargout = {Fp,Fm,F0};
        end
    end

end