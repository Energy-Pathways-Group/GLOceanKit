classdef TestRandomFlow < matlab.unittest.TestCase
    properties
        wvt
        solutionGroup
    end

    properties (ClassSetupParameter)
        % Lxyz = struct('Lxyz',[15e3, 15e3, 1300]);
        % Nxyz = struct('Nx8Ny8Nz5',[8 8 5]);
        % Nxyz = struct('Nx16Ny16Nz5',[16 16 5]);
        % transform = {'constant','hydrostatic','boussinesq'};
        Lxyz = struct('Lxyz',[1000, 500, 500]);
        Nxyz = struct('Nx16Ny8Nz9',[16 8 9]);
        % Nxyz = struct('Nx32N16Nz17',[32 16 17]);
        transform = {'constant'};
        % orthogonalSolutionGroup = {'WVInertialOscillationComponent','WVMeanDensityAnomalyComponent','WVInternalGravityWaveComponent','WVGeostrophicComponent'}
        % flowComponent = {'WVMeanDensityAnomalyComponent'}
    end

    methods (TestClassSetup)
        function classSetup(testCase,Lxyz,Nxyz,transform)
            switch transform
                case 'constant'
                    testCase.wvt = WVTransformConstantStratification(Lxyz, Nxyz, latitude=33, isHydrostatic=0);
                case 'hydrostatic'
                    testCase.wvt = WVTransformHydrostatic(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
                case 'boussinesq'
                    testCase.wvt = WVTransformBoussinesq(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
            end
        end
    end

    methods (TestParameterDefinition,Static)
        function [flowComponent] = initializeProperty(Lxyz,Nxyz,transform)
            % If you want to dynamically adjust the test parameters, you
            % have to do it here.
            switch transform
                case 'constant'
                    tmpwvt = WVTransformConstantStratification(Lxyz, Nxyz, latitude=33, isHydrostatic=0);
                case 'hydrostatic'
                    tmpwvt = WVTransformHydrostatic(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
                case 'boussinesq'
                    tmpwvt = WVTransformBoussinesq(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
            end
            
            flowComponent = tmpwvt.flowComponentNames;
        end
    end

    properties (TestParameter)
        flowComponent
    end

    methods (Test)
        function testSolution(self,flowComponent)
            self.wvt.initWithRandomFlow(flowComponent,uvMax=0.1);
            self.verifyEqual(self.wvt.totalEnergy,self.wvt.totalEnergySpatiallyIntegrated, "AbsTol",1e-7,"RelTol",1e-7);
        end

    end

end