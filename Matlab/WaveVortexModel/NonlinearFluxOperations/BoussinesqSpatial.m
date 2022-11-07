classdef BoussinesqSpatial < WVNonlinearFluxOperation
    properties
        dLnN2
    end
    methods
        function self = BoussinesqSpatial(wvt)
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            fluxVar(1) = WVVariableAnnotation('Fp',{'k','l','j'},'m/s2', 'non-linear flux into Ap');
            fluxVar(2) = WVVariableAnnotation('Fm',{'k','l','j'},'m/s2', 'non-linear flux into Am');
            fluxVar(3) = WVVariableAnnotation('F0',{'k','l','j'},'m/s', 'non-linear flux into A0');

            self@WVNonlinearFluxOperation('BoussinesqSpatial',fluxVar);

            if isa(wvt,'WVTransformConstantStratification')
                self.dLnN2 = zeros(size(wvt.z));
            elseif isa(wvt,'WVTransformHydrostatic')
                self.dLnN2 = wvt.dLnN2;
            end
        end

        function varargout = compute(self,wvt,varargin)
            varargout = cell(1,self.nVarOut);
            uNL = wvt.u .* wvt.diffX(wvt.u)   + wvt.v .* wvt.diffY(wvt.u)   + wvt.w .*  wvt.diffZF(wvt.u);
            vNL = wvt.u .* wvt.diffX(wvt.v)   + wvt.v .* wvt.diffY(wvt.v)   + wvt.w .*  wvt.diffZF(wvt.v);
            nNL = wvt.u .* wvt.diffX(wvt.eta) + wvt.v .* wvt.diffY(wvt.eta) + wvt.w .* (wvt.diffZG(wvt.eta) + wvt.eta .* shiftdim(self.dLnN2,-2));

            [varargout{:}] = wvt.transformUVEtaToWaveVortex(uNL,vNL,nNL,wvt.t);
        end
    end

end