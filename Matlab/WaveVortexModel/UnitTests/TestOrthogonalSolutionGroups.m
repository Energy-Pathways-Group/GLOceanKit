classdef TestOrthogonalSolutionGroups < matlab.unittest.TestCase
    properties
        wvt
        solutionGroup
    end

    properties (ClassSetupParameter)
        % transform = {'constant','hydrostatic','boussinesq'};
        transform = {'boussinesq'};
        orthogonalSolutionGroup = {'WVInertialOscillationSolutionGroup'}
    end

    methods (TestClassSetup)
        function classSetup(testCase,transform,orthogonalSolutionGroup)
            switch transform
                case 'constant'
                    testCase.wvt = WVTransformConstantStratification([1, 10, 4], [16 16 9]);
                case 'hydrostatic'
                    testCase.wvt = WVTransformHydrostatic([1, 10, 4], [16 16 9], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
                case 'boussinesq'
                    testCase.wvt = WVTransformBoussinesq([1, 10, 4], [16 16 9], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
            end
            switch orthogonalSolutionGroup
                case 'WVInertialOscillationSolutionGroup'
                    testCase.solutionGroup = WVInertialOscillationSolutionGroup(testCase.wvt);
            end
        end
    end

    methods (TestParameterDefinition,Static)
        function [solutionIndex] = initializeProperty(transform,orthogonalSolutionGroup)
            % If you want to dynamically adjust the test parameters, you
            % have to do it here.
            switch transform
                case 'constant'
                    tmpwvt = WVTransformConstantStratification([1, 10, 4], [16 16 9]);
                case 'hydrostatic'
                    tmpwvt = WVTransformHydrostatic([1, 10, 4], [16 16 9], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
                case 'boussinesq'
                    tmpwvt = WVTransformBoussinesq([1, 10, 4], [16 16 9], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
            end
            switch orthogonalSolutionGroup
                case 'WVInertialOscillationSolutionGroup'
                    solnGroup = WVInertialOscillationSolutionGroup(tmpwvt);
            end
            for iSoln = 1:solnGroup.nUniqueSolutions
                solutionIndex.(sprintf('solution_%d',iSoln)) = iSoln;
            end
        end
    end

    properties (TestParameter)
        solutionIndex
    end

    methods (Test)
        function testSolution(self,solutionIndex)
            wvt = self.wvt;
            soln = self.solutionGroup.uniqueSolutionAtIndex(solutionIndex,amplitude='random');
            wvt.initWithUVEta(soln.u(wvt.X,wvt.Y,wvt.Z,wvt.t), soln.v(wvt.X,wvt.Y,wvt.Z,wvt.t),soln.eta(wvt.X,wvt.Y,wvt.Z,wvt.t));

            self.verifyEqual(wvt.u,soln.u(wvt.X,wvt.Y,wvt.Z,wvt.t), "AbsTol",1e-7,"RelTol",1e-7);

        end


    end

end