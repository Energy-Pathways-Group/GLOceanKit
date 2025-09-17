classdef TestRadialTransformation < matlab.unittest.TestCase
    % This test will create randomized flow components, compute total
    % energy and enstrophy variance in the full domain and radial
    % wavenumber domain
    properties
        wvt
    end

    properties (ClassSetupParameter)
        Lxyz = struct('Lxyz',[15e3 15e3 1300]);
        Nxyz = struct('Nx32Ny32Nz17',[32 32 17]);
        transform = {'hydrostatic'};
    end

    methods (TestClassSetup)
        function classSetup(testCase,Lxyz,Nxyz,transform)
            switch transform
                case 'constant'
                    testCase.wvt = WVTransformConstantStratification(Lxyz, Nxyz);
                case 'hydrostatic'
                    testCase.wvt = WVTransformHydrostatic(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
                case 'boussinesq'
                    testCase.wvt = WVTransformBoussinesq(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
            end
        end
    end

    properties (TestParameter)
        flowComponent = {'geostrophic','mda','wave','inertial'}
    end

    methods (Test)
        function testRadialWavenumberVariance(self,flowComponent)
            self.wvt.initWithRandomFlow(flowComponent);
            varianceMatrix = abs(self.wvt.Ap).^2 + abs(self.wvt.Am).^2 + abs(self.wvt.A0).^2;
            radialVarianceMatrix = self.wvt.transformToRadialWavenumber(varianceMatrix);

            self.verifyEqual(sum(radialVarianceMatrix(:)),sum(radialVarianceMatrix(:)), "AbsTol",1e-7,"RelTol",1e-7);
        end


    end

end