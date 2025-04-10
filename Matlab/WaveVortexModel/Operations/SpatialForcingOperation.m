classdef SpatialForcingOperation < WVOperation
    % Computes the nonlinear flux for a WVTransform

    properties
        Fpv
    end

    methods

        function self = SpatialForcingOperation(wvt)
            arguments
                wvt WVTransform
            end
            outputVariables = WVVariableAnnotation.empty(0,0);
            for i=1:length(wvt.forcing)

                if isa(wvt,"WVTransformBarotropicQG") || isa(wvt,"WVTransformStratifiedQG")
                    name = replace(replace(join( ["Fqgpv_", string(wvt.forcing(i).name)],"")," ","_"),"-","_");
                    outputVariables(i) = WVVariableAnnotation(name,wvt.spatialDimensionNames(),'s^{-2}', join(['spatial representation of qgpv forcing',string(wvt.forcing(i).name)]));
                elseif isa(wvt,"WVTransformHydrostatic")
                    name = replace(replace(join( ["Fu_", string(wvt.forcing(i).name)],"")," ","_"),"-","_");
                    outputVariables((i-1)*3+1) = WVVariableAnnotation(name,wvt.spatialDimensionNames(),'m s^{-2}', join(['spatial representation of hydrostatic forcing on the x-momentum equation',string(wvt.forcing(i).name)]));

                    name = replace(replace(join( ["Fv_", string(wvt.forcing(i).name)],"")," ","_"),"-","_");
                    outputVariables((i-1)*3+2) = WVVariableAnnotation(name,wvt.spatialDimensionNames(),'m s^{-2}', join(['spatial representation of hydrostatic forcing on the y-momentum equation',string(wvt.forcing(i).name)]));

                    name = replace(replace(join( ["Feta_", string(wvt.forcing(i).name)],"")," ","_"),"-","_");
                    outputVariables((i-1)*3+3) = WVVariableAnnotation(name,wvt.spatialDimensionNames(),'m s^{-1}', join(['spatial representation of hydrostatic forcing on the scaled density perturbation equation',string(wvt.forcing(i).name)]));
                end
            end
            self@WVOperation('spatial forcing',outputVariables,@disp);
            self.Fpv = zeros(wvt.spatialMatrixSize);
        end

        function varargout = compute(self,wvt,varargin)
            varargout = cell(self.nVarOut,1);

            if isa(wvt,"WVTransformBarotropicQG") || isa(wvt,"WVTransformStratifiedQG")
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
            elseif isa(wvt,"WVTransformHydrostatic")
                Fu=0;Fv=0;Feta=0; % this isn't good, need to cached
                for i=1:length(self.spatialFluxForcing)
                    Fu0=Fu;Fv0=Fv;Feta0=Feta;
                    [Fu, Fv, Feta] = self.spatialFluxForcing(i).addHydrostaticSpatialForcing(self, Fu, Fv, Feta);
                    [Fp,Fm,F0] = self.transformUVEtaToWaveVortex(Fu-Fu0, Fv-Fv0, Feta-Feta0,self.t);
                    F{self.spatialFluxForcing(i).name} = struct("Fp",Fp,"Fm",Fm,"F0",F0);
                end

            end
        end

    end

end