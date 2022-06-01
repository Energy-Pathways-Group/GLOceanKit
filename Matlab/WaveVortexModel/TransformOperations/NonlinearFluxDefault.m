classdef NonlinearFluxDefault < TransformOperation

    methods
        function self = NonlinearFluxDefault()
            fluxVar(1) = StateVariable('Fp',{'k','l','j'},'m/s2', 'non-linear flux into Ap');
            fluxVar(2) = StateVariable('Fm',{'k','l','j'},'m/s2', 'non-linear flux into Am');
            fluxVar(3) = StateVariable('F0',{'k','l','j'},'m/s', 'non-linear flux into A0');

            self@TransformOperation('NonlinearFluxDefault',fluxVar);
        end

        function varargout = Compute(self,wvt,varargin)
            phase = exp(self.iOmega*(wvt.t-wvt.t0));
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
            [U,Ux,Uy,Uz] = wvt.TransformToSpatialDomainWithFAllDerivatives(Ubar);
            [V,Vx,Vy,Vz] = wvt.TransformToSpatialDomainWithFAllDerivatives(Vbar);
            W = wvt.TransformToSpatialDomainWithG(Wbar);
            [ETA,ETAx,ETAy,ETAz] = wvt.TransformToSpatialDomainWithGAllDerivatives(Nbar);

            % Compute the nonlinear terms in the spatial domain
            % (pseudospectral!)
            uNL = -U.*Ux - V.*Uy - W.*Uz;
            vNL = -U.*Vx - V.*Vy - W.*Vz;
            nNL = -U.*ETAx - V.*ETAy - W.*(ETAz + ETA.*shiftdim(self.dLnN2,-2));

            % Now apply the operator S^{-1} and then T_\omega^{-1}
            uNLbar = wvt.TransformFromSpatialDomainWithF(uNL);
            vNLbar = wvt.TransformFromSpatialDomainWithF(vNL);
            nNLbar = wvt.TransformFromSpatialDomainWithG(nNL);

            Fp = (self.ApU.*uNLbar + self.ApV.*vNLbar + self.ApN.*nNLbar) .* conj(phase);
            Fm = (self.AmU.*uNLbar + self.AmV.*vNLbar + self.AmN.*nNLbar) .* phase;
            F0 = (self.A0U.*uNLbar + self.A0V.*vNLbar + self.A0N.*nNLbar);

            varargout = {Fp,Fm,F0};
        end
    end

end