classdef TestRadialTransformation < matlab.unittest.TestCase
    % This test will create randomized flow components, compute total
    % energy and enstrophy variance in the full domain and radial
    % wavenumber domain
    properties
        wvt
        flowComponent
    end

    properties (ClassSetupParameter)
        Lxyz = struct('Lxyz',[15e3 15e3 1300]);
        Nxyz = struct('Nx32Ny32Nz17',[32 32 17]);
        transform = {'constant'};
        flowComponentType = {'WVInertialOscillationComponent','WVMeanDensityAnomalyComponent','WVInternalGravityWaveComponent','WVGeostrophicComponent'}
    end

    methods (TestClassSetup)
        function classSetup(testCase,Lxyz,Nxyz,transform,flowComponentType)
            switch transform
                case 'constant'
                    testCase.wvt = WVTransformConstantStratification(Lxyz, Nxyz);
                case 'hydrostatic'
                    testCase.wvt = WVTransformHydrostatic(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
                case 'boussinesq'
                    testCase.wvt = WVTransformBoussinesq(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
            end
            switch flowComponentType
                case 'WVInertialOscillationComponent'
                    testCase.flowComponent = WVInertialOscillationComponent(testCase.wvt);
                case 'WVMeanDensityAnomalyComponent'
                    testCase.flowComponent = WVMeanDensityAnomalyComponent(testCase.wvt);
                case 'WVInternalGravityWaveComponent'
                    testCase.flowComponent = WVInternalGravityWaveComponent(testCase.wvt);
                case 'WVGeostrophicComponent'
                    testCase.flowComponent = WVGeostrophicComponent(testCase.wvt);
            end
        end
    end

    methods (Test)
        function testRadialWavenumberVariance(self)
            varianceMatrix = abs(self.wvt.Ap).^2 + abs(self.wvt.Am).^2 + abs(self.wvt.A0).^2;
            radialVarianceMatrix = self.wvt.transformToRadialWavenumber(varianceMatrix);

            self.verifyEqual(sum(radialVarianceMatrix(:)),sum(radialVarianceMatrix(:)), "AbsTol",1e-7,"RelTol",1e-7);
        end


    end

end