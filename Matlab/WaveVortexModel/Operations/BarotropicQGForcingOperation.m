classdef BarotropicQGForcingOperation < WVOperation
    % Computes the nonlinear flux for a WVTransform

    properties
        Fpv
    end

    methods

        function self = BarotropicQGForcingOperation(wvt)
            arguments
                wvt WVTransform
            end
            outputVariables = WVVariableAnnotation.empty(0,0);
            for i=1:length(wvt.forcing)
                name = join( ["force_", string(wvt.forcing(i).name)],"");
                name = replace(name," ","_");
                name = replace(name,"-","_");
                outputVariables(i) = WVVariableAnnotation(name,WVTransformBarotropicQG.spatialDimensionNames(),'s^{-2}', join(['spatial representation of qgpv forcing',string(wvt.forcing(i).name)]));
            end
            self@WVOperation('spatial forcing',outputVariables,@disp);
            self.Fpv = zeros(wvt.spatialMatrixSize);
        end

        function varargout = compute(self,wvt,varargin)
            varargout = cell(self.nVarOut,1);
            self.Fpv = 0*self.Fpv;
            iForce = 0;
            for i=1:length(wvt.spatialFluxForcing)
                iForce = iForce + 1;
                Fpv0 = self.Fpv;
                self.Fpv = wvt.spatialFluxForcing(i).addPotentialVorticitySpatialForcing(wvt,self.Fpv);
                varargout{iForce} = (self.Fpv-Fpv0);
            end
            F0 = wvt.transformFromSpatialDomainWithFourier(self.Fpv);
            for i=1:length(wvt.spectralFluxForcing)
               iForce = iForce + 1;
               F0_i = F0;
               F0 = wvt.spectralFluxForcing(i).addPotentialVorticitySpectralForcing(wvt,F0);
               varargout{iForce} = wvt.transformToSpatialDomainWithFourier(F0 - F0_i);
            end
            for i=1:length(wvt.spectralAmplitudeForcing)
               iForce = iForce + 1;
               F0_i = F0;
               F0 = wvt.spectralAmplitudeForcing(i).setPotentialVorticitySpectralForcing(wvt,F0);
               varargout{iForce} = wvt.transformToSpatialDomainWithFourier(F0 - F0_i);
            end
        end

    end

end