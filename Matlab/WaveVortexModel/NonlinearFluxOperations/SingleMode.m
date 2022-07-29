classdef SingleMode < WVNonlinearFluxOperation

    methods
        function self = SingleMode()
            fluxVar(1) = WVVariableAnnotation('Fp',{'k','l','j'},'m/s2', 'non-linear flux into Ap');
            fluxVar(2) = WVVariableAnnotation('Fm',{'k','l','j'},'m/s2', 'non-linear flux into Am');
            fluxVar(3) = WVVariableAnnotation('F0',{'k','l','j'},'m/s', 'non-linear flux into A0');

            self@WVNonlinearFluxOperation('SingleMode',fluxVar);
        end

        function varargout = compute(self,wvt,varargin)
            varargout = cell(1,self.nVarOut);
            uNL = wvt.u .* wvt.diffX(wvt.u)   + wvt.v .* wvt.diffY(wvt.u)  ;
            vNL = wvt.u .* wvt.diffX(wvt.v)   + wvt.v .* wvt.diffY(wvt.v)  ;
            nNL = wvt.u .* wvt.diffX(wvt.eta) + wvt.v .* wvt.diffY(wvt.eta);

            [varargout{:}] = wvt.transformUVEtaToWaveVortex(uNL,vNL,nNL,wvt.t);
        end
    end

end