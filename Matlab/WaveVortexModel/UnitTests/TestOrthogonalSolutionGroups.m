classdef TestOrthogonalSolutionGroups < matlab.unittest.TestCase
    properties
        wvt
        solutionGroup
    end

    properties (ClassSetupParameter)
        Lxyz = struct('Lxyz',[15e3, 15e3, 1300]);
        Nxyz = struct('Nx8Ny8Nz5',[8 8 5]);
        % Nxyz = struct('Nx16Ny16Nz5',[16 16 5]);
        transform = {'constant','hydrostatic','boussinesq'};
        % transform = {'hydrostatic'};
        orthogonalSolutionGroup = {'WVInertialOscillationComponent','WVMeanDensityAnomalyComponent','WVInternalGravityWaveComponent','WVGeostrophicComponent'}
        % orthogonalSolutionGroup = {'WVInternalGravityWaveComponent'}
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
                case 'WVInertialOscillationComponent'
                    testCase.solutionGroup = WVInertialOscillationComponent(testCase.wvt);
                case 'WVMeanDensityAnomalyComponent'
                    testCase.solutionGroup = WVMeanDensityAnomalyComponent(testCase.wvt);
                case 'WVInternalGravityWaveComponent'
                    testCase.solutionGroup = WVInternalGravityWaveComponent(testCase.wvt);
                case 'WVGeostrophicComponent'
                    testCase.solutionGroup = WVGeostrophicComponent(testCase.wvt);
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
                case 'WVInertialOscillationComponent'   
                    solnGroup = WVInertialOscillationComponent(tmpwvt);
                case 'WVMeanDensityAnomalyComponent'
                    solnGroup = WVMeanDensityAnomalyComponent(tmpwvt);
                case 'WVInternalGravityWaveComponent'
                    solnGroup = WVInternalGravityWaveComponent(tmpwvt);
                case 'WVGeostrophicComponent'
                    solnGroup = WVGeostrophicComponent(tmpwvt);
            end
            for iSoln = 1:solnGroup.nModes
                solutionIndex.(sprintf('solution_%d',iSoln)) = iSoln;
            end
        end
    end

    properties (TestParameter)
        solutionIndex
    end

    methods (Test)
        function testSolution(self,solutionIndex)
            self.wvt.t=0;
            args = {self.wvt.X,self.wvt.Y,self.wvt.Z,self.wvt.t};
            soln = self.solutionGroup.solutionForModeAtIndex(solutionIndex,amplitude='random');
            self.wvt.initWithUVEta(soln.u(args{:}), soln.v(args{:}),soln.eta(args{:}));

            % Advance forward in time to confirm that phases are correctly
            % changing.
            self.wvt.t = 86400;
            args = {self.wvt.X,self.wvt.Y,self.wvt.Z,self.wvt.t};
            self.verifyThat(self.wvt.u,IsSameSolutionAs(soln.u(args{:})),'u');
            self.verifyThat(self.wvt.v,IsSameSolutionAs(soln.v(args{:})),'v');
            self.verifyThat(self.wvt.w,IsSameSolutionAs(soln.w(args{:})),'w');
            self.verifyThat(self.wvt.eta,IsSameSolutionAs(soln.eta(args{:})),'eta');
            self.verifyThat(self.wvt.p,IsSameSolutionAs(soln.p(args{:})),'p');
            self.verifyThat(self.wvt.qgpv,IsSameSolutionAs(soln.qgpv(args{:})),'qgpv');

            self.verifyEqual(self.wvt.totalEnergy,soln.depthIntegratedTotalEnergy(isHydrostatic=self.wvt.isHydrostatic), "AbsTol",1e-7,"RelTol",1e-3);
            self.verifyEqual(self.wvt.totalEnstrophy,soln.depthIntegratedTotalEnstrophy, "AbsTol",1e-7,"RelTol",1e-3);
        end

        function testFTransformAllDerivatives(self,solutionIndex)
            args = {self.wvt.X,self.wvt.Y,self.wvt.Z,self.wvt.t};
            soln = self.solutionGroup.solutionForModeAtIndex(solutionIndex,amplitude='random');
            self.wvt.initWithUVEta(soln.u(args{:}), soln.v(args{:}),soln.eta(args{:}));
            [U,Ux,Uy,Uz] = self.wvt.transformToSpatialDomainWithFAllDerivatives(Apm=self.wvt.UAp.*self.wvt.Apt + self.wvt.UAm.*self.wvt.Amt,A0=self.wvt.UA0.*self.wvt.A0t);

            self.verifyThat(U,IsSameSolutionAs(soln.u(args{:})),'FAllDerivatives-U');
            self.verifyThat(Ux,IsSameSolutionAs(self.wvt.diffX(self.wvt.u)),'FAllDerivatives-Ux');
            self.verifyThat(Uy,IsSameSolutionAs(self.wvt.diffY(self.wvt.u)),'FAllDerivatives-Uy');
            self.verifyThat(Uz,IsSameSolutionAs(self.wvt.diffZF(self.wvt.u)),'FAllDerivatives-Uz');
        end

        function testGTransformAllDerivatives(self,solutionIndex)
            args = {self.wvt.X,self.wvt.Y,self.wvt.Z,self.wvt.t};
            soln = self.solutionGroup.solutionForModeAtIndex(solutionIndex,amplitude='random');
            self.wvt.initWithUVEta(soln.u(args{:}), soln.v(args{:}),soln.eta(args{:}));
            [ETA,ETAx,ETAy,ETAz] = self.wvt.transformToSpatialDomainWithGAllDerivatives(Apm=self.wvt.NAp.*self.wvt.Apt + self.wvt.NAm.*self.wvt.Amt,A0=self.wvt.NA0.*self.wvt.A0t);

            self.verifyThat(ETA,IsSameSolutionAs(soln.eta(args{:})),'GAllDerivatives-ETA');
            self.verifyThat(ETAx,IsSameSolutionAs(self.wvt.diffX(self.wvt.eta)),'GAllDerivatives-ETAx');
            self.verifyThat(ETAy,IsSameSolutionAs(self.wvt.diffY(self.wvt.eta)),'GAllDerivatives-ETAy');
            self.verifyThat(ETAz,IsSameSolutionAs(self.wvt.diffZG(self.wvt.eta)),'GAllDerivatives-ETAz');
        end

    end

end