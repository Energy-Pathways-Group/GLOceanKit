classdef BoussinesqConstantN < WVNonlinearFluxOperation

    properties
        shouldAntialias = 0
        AA
        shouldDamp = 0
        nu_xy
        nu_z
    end

    methods
        function self = BoussinesqConstantN(wvt,options)
            arguments
                wvt WVTransform {mustBeNonempty}
                options.uv_damp (1,1) double % characteristic speed used to set the damping. Try using uMax.
                options.w_damp (1,1) double % characteristic speed used to set the damping. Try using wMax
                options.nu (1,1) double
                options.shouldAntialias double = 0
            end
            fluxVar(1) = WVVariableAnnotation('Fp',{'k','l','j'},'m/s2', 'non-linear flux into Ap',detailedDescription='- topic: State Variables');
            fluxVar(2) = WVVariableAnnotation('Fm',{'k','l','j'},'m/s2', 'non-linear flux into Am',detailedDescription='- topic: State Variables');
            fluxVar(3) = WVVariableAnnotation('F0',{'k','l','j'},'m/s', 'non-linear flux into A0',detailedDescription='- topic: State Variables');

            self@WVNonlinearFluxOperation('BoussinesqConstantN',fluxVar);
            self.shouldAntialias = options.shouldAntialias;
            
            if self.shouldAntialias == 1
                self.AA = ~(wvt.maskForAliasedModes(jFraction=2/3));
            end
            
            if isfield(options,'u_damp')
                if self.shouldAntialias == 1
                    self.nu_xy = (3/2)*(wvt.x(2)-wvt.x(1))*options.uv_damp/(pi^2);
                else
                    self.nu_xy = (wvt.x(2)-wvt.x(1))*options.uv_damp/(pi^2);
                end
            end

            if isfield(options,'w_damp')
                if self.shouldAntialias == 1
                    self.nu_z = (3/2)*(wvt.z(2)-wvt.z(1))*options.w_damp/(pi^2);
                else
                    self.nu_z = (wvt.z(2)-wvt.z(1))*options.w_damp/(pi^2);
                end
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
            [~,ETAx,ETAy,ETAz] = wvt.transformToSpatialDomainWithGAllDerivatives(Nbar);

            % compute the nonlinear terms in the spatial domain
            % (pseudospectral!)
            uNL = -U.*Ux - V.*Uy - W.*Uz;
            vNL = -U.*Vx - V.*Vy - W.*Vz;
            nNL = -U.*ETAx - V.*ETAy - W.*ETAz;

            % Now apply the operator S^{-1} and then T_\omega^{-1}
            uNLbar = wvt.transformFromSpatialDomainWithF(uNL);
            vNLbar = wvt.transformFromSpatialDomainWithF(vNL);
            nNLbar = wvt.transformFromSpatialDomainWithG(nNL);

            Fp = (wvt.ApU.*uNLbar + wvt.ApV.*vNLbar + wvt.ApN.*nNLbar) .* conj(phase);
            Fm = (wvt.AmU.*uNLbar + wvt.AmV.*vNLbar + wvt.AmN.*nNLbar) .* phase;
            F0 = (wvt.A0U.*uNLbar + wvt.A0V.*vNLbar + wvt.A0N.*nNLbar);

            varargout = {Fp,Fm,F0};
        end
    end

end