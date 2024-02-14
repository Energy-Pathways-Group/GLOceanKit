classdef TestOrthogonalSolutionGroups < matlab.unittest.TestCase
    properties
        wvt
        solutionGroup
    end

    properties (ClassSetupParameter)
        Lxyz = struct('Lxyz',[1 10 4]);
        Nxyz = struct('Nx8Ny8Nz5',[8 8 5]);
        %transform = {'constant','hydrostatic','boussinesq'};
        transform = {'hydrostatic'};
        orthogonalSolutionGroup = {'WVInertialOscillationSolutionGroup','WVMeanDensityAnomalySolutionGroup','WVInternalGravityWaveSolutionGroup','WVGeostrophicSolutionGroup'}
    end

    methods (TestClassSetup)
        function classSetup(testCase,Lxyz,Nxyz,transform,orthogonalSolutionGroup)
            switch transform
                case 'constant'
                    testCase.wvt = WVTransformConstantStratification(Lxyz, Nxyz);
                case 'hydrostatic'
                    testCase.wvt = WVTransformHydrostatic(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
                case 'boussinesq'
                    testCase.wvt = WVTransformBoussinesq(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
            end
            switch orthogonalSolutionGroup
                case 'WVInertialOscillationSolutionGroup'
                    testCase.solutionGroup = WVInertialOscillationSolutionGroup(testCase.wvt);
                case 'WVMeanDensityAnomalySolutionGroup'
                    testCase.solutionGroup = WVMeanDensityAnomalySolutionGroup(testCase.wvt);
                case 'WVInternalGravityWaveSolutionGroup'
                    testCase.solutionGroup = WVInternalGravityWaveSolutionGroup(testCase.wvt);
                case 'WVGeostrophicSolutionGroup'
                    testCase.solutionGroup = WVGeostrophicSolutionGroup(testCase.wvt);
            end
        end
    end

    methods (TestParameterDefinition,Static)
        function [solutionIndex] = initializeProperty(Lxyz,Nxyz,transform,orthogonalSolutionGroup)
            % If you want to dynamically adjust the test parameters, you
            % have to do it here.
            switch transform
                case 'constant'
                    tmpwvt = WVTransformConstantStratification(Lxyz, Nxyz);
                case 'hydrostatic'
                    tmpwvt = WVTransformHydrostatic(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
                case 'boussinesq'
                    tmpwvt = WVTransformBoussinesq(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
            end
            % ideally we re-write the initialization to only depend
            % on Lxyz and Nyxz
            switch orthogonalSolutionGroup
                case 'WVInertialOscillationSolutionGroup'   
                    solnGroup = WVInertialOscillationSolutionGroup(tmpwvt);
                case 'WVMeanDensityAnomalySolutionGroup'
                    solnGroup = WVMeanDensityAnomalySolutionGroup(tmpwvt);
                case 'WVInternalGravityWaveSolutionGroup'
                    solnGroup = WVInternalGravityWaveSolutionGroup(tmpwvt);
                case 'WVGeostrophicSolutionGroup'
                    solnGroup = WVGeostrophicSolutionGroup(tmpwvt);
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
            args = {self.wvt.X,self.wvt.Y,self.wvt.Z,self.wvt.t};
            soln = self.solutionGroup.uniqueSolutionAtIndex(solutionIndex,amplitude='random');
            self.wvt.initWithUVEta(soln.u(args{:}), soln.v(args{:}),soln.eta(args{:}));

            self.verifyThat(self.wvt.u,IsSameSolutionAs(soln.u(args{:})),'u');
            self.verifyThat(self.wvt.v,IsSameSolutionAs(soln.v(args{:})),'v');
            self.verifyThat(self.wvt.w,IsSameSolutionAs(soln.w(args{:})),'w');
            self.verifyThat(self.wvt.eta,IsSameSolutionAs(soln.eta(args{:})),'eta');
            self.verifyThat(self.wvt.p,IsSameSolutionAs(soln.p(args{:})),'p');
        end

        function testFTransformAllDerivatives(self,solutionIndex)
            args = {self.wvt.X,self.wvt.Y,self.wvt.Z,self.wvt.t};
            soln = self.solutionGroup.uniqueSolutionAtIndex(solutionIndex,amplitude='random');
            self.wvt.initWithUVEta(soln.u(args{:}), soln.v(args{:}),soln.eta(args{:}));
            [U,Ux,Uy,Uz] = self.wvt.transformToSpatialDomainWithFAllDerivatives(Apm=self.wvt.UAp.*self.wvt.Apt + self.wvt.UAm.*self.wvt.Amt,A0=self.wvt.UA0.*self.wvt.A0t);

            self.verifyThat(U,IsSameSolutionAs(soln.u(args{:})),'FAllDerivatives-U');
            self.verifyThat(Ux,IsSameSolutionAs(self.wvt.diffX(self.wvt.u)),'FAllDerivatives-Ux');
            self.verifyThat(Uy,IsSameSolutionAs(self.wvt.diffY(self.wvt.u)),'FAllDerivatives-Uy');
            self.verifyThat(Uz,IsSameSolutionAs(self.wvt.diffZF(self.wvt.u)),'FAllDerivatives-Uz');
        end

        function testGTransformAllDerivatives(self,solutionIndex)
            args = {self.wvt.X,self.wvt.Y,self.wvt.Z,self.wvt.t};
            soln = self.solutionGroup.uniqueSolutionAtIndex(solutionIndex,amplitude='random');
            self.wvt.initWithUVEta(soln.u(args{:}), soln.v(args{:}),soln.eta(args{:}));
            [ETA,ETAx,ETAy,ETAz] = self.wvt.transformToSpatialDomainWithGAllDerivatives(Apm=self.wvt.NAp.*self.wvt.Apt + self.wvt.NAm.*self.wvt.Amt,A0=self.wvt.NA0.*self.wvt.A0t);

            self.verifyThat(ETA,IsSameSolutionAs(soln.eta(args{:})),'GAllDerivatives-ETA');
            self.verifyThat(ETAx,IsSameSolutionAs(self.wvt.diffX(self.wvt.eta)),'GAllDerivatives-ETAx');
            self.verifyThat(ETAy,IsSameSolutionAs(self.wvt.diffY(self.wvt.eta)),'GAllDerivatives-ETAy');
            self.verifyThat(ETAz,IsSameSolutionAs(self.wvt.diffZG(self.wvt.eta)),'GAllDerivatives-ETAz');
        end

    end

end