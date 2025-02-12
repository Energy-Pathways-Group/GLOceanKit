classdef TestInternalGravityWaveMethods < matlab.unittest.TestCase
    properties
        wvt
        solutionGroup
    end

    properties (ClassSetupParameter)
        % Lxyz = struct('Lxyz',[750e0, 750e0, 375]);
        Lxyz = struct('Lxyz',[750e2, 750e2, 1300]);
        %Nxyz = struct('Nx64Ny64Nz30',[64 64 32]);
        % Nxyz = struct('Nx16Ny16Nz5',[16 16 5]);
        Nxyz = struct('Nx32Ny32Nz10',2*[16 16 5]);
        %transform = {'constant','hydrostatic','boussinesq'};
        transform = {'hydrostatic'};
    end

    methods (TestClassSetup)
        function classSetup(testCase,Lxyz,Nxyz,transform)
            switch transform
                case 'constant'
                    testCase.wvt = WVTransformConstantStratification(Lxyz, Nxyz);
                case 'hydrostatic'
                    testCase.wvt = WVTransformHydrostatic(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)), shouldAntialias=1);
                case 'boussinesq'
                    testCase.wvt = WVTransformBoussinesq(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)), shouldAntialias=1);
            end
            testCase.solutionGroup = WVGeostrophicComponent(testCase.wvt);
        end
    end

    methods (Test)
        function testRemoveAllWaves(self)
            % In this test we intialize with a random flow state, confirm
            % that both total energy and geostrophic energy are present,
            % remove all the geostrophic energy, and then confirm that there
            % is no geostrophic energy remaining, and that the total energy is
            % the same as the initial total, minus the initial geostrophic.
            self.wvt.removeAll();

            % currently we need to make this small-ish to avoid density
            % overturns.
            self.wvt.initWithRandomFlow(uvMax=0.01);

            initialTotalEnergy = self.wvt.totalEnergy;
            initialWaveEnergy = self.wvt.waveEnergy;
            self.verifyGreaterThan(initialTotalEnergy,0.0);
            self.verifyGreaterThan(initialWaveEnergy,0.0);

            self.wvt.removeAllWaves();
            finalTotalEnergy = self.wvt.totalEnergy;
            finalWaveEnergy = self.wvt.waveEnergy;

            self.verifyEqual(finalWaveEnergy,0.0);
            self.verifyEqual(finalTotalEnergy,initialTotalEnergy-initialWaveEnergy, "AbsTol",1e-7,"RelTol",1e-7);
        end

        % function testDensityBoundsViolation(self)
        %     % confirm an error when the amplitude is too large
        %     self.wvt.removeAll();
        % 
        %     Le = 120e3;
        %     He = 200;
        %     x0 = (1/2)*max(self.wvt.x); y0=max(self.wvt.y)/2;
        % 
        %     U = 0.20; % m/s
        %     psibar = @(z) (pi*Le*Le/(self.wvt.Lx*self.wvt.Ly))*U*(Le/sqrt(2))*exp(1/2)*exp(-(z/He/sqrt(2)).^2 );
        %     psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2 -(z/He/sqrt(2)).^2 ) - psibar(z);
        %     self.verifyError(@() self.wvt.initWithGeostrophicStreamfunction(psi),'WVTransform:DensityBoundsViolation' );
        % 
        %     self.wvt.removeAll();
        %     U = -1.00; % m/s
        %     psibar = @(z) (pi*Le*Le/(self.wvt.Lx*self.wvt.Ly))*U*(Le/sqrt(2))*exp(1/2)*exp(-(z/He/sqrt(2)).^2 );
        %     psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2 -(z/He/sqrt(2)).^2 ) - psibar(z);
        %     self.verifyError(@() self.wvt.initWithGeostrophicStreamfunction(psi),'WVTransform:DensityBoundsViolation' );
        % end

        function testSetWaveModes(self)
            self.wvt.removeAll();

            kMode = 3; lMode = 5; j=3; phi = pi*0.3; u=0.05; sign = +1;

            self.wvt.t = 22654;
            self.wvt.setWaveModes(kMode=kMode,lMode=lMode,j=j,phi=phi,u=u,sign=sign);

            soln = self.wvt.waveComponent.internalGravityWaveSolution(kMode,lMode,j,u,phi,sign,amplitudeIsMaxU=1,t=self.wvt.t);

            self.wvt.t = 86400;
            args = {self.wvt.X,self.wvt.Y,self.wvt.Z,self.wvt.t};
            self.verifyThat(self.wvt.u,IsSameSolutionAs(soln.u(args{:})),'u');
            self.verifyThat(self.wvt.v,IsSameSolutionAs(soln.v(args{:})),'v');
            self.verifyThat(self.wvt.w,IsSameSolutionAs(soln.w(args{:})),'w');
            self.verifyThat(self.wvt.eta,IsSameSolutionAs(soln.eta(args{:})),'eta');
            self.verifyThat(self.wvt.p,IsSameSolutionAs(soln.p(args{:})),'p');
            self.verifyThat(self.wvt.qgpv,IsSameSolutionAs(soln.qgpv(args{:})),'qgpv');

            self.verifyEqual(self.wvt.totalEnergy,soln.depthIntegratedTotalEnergy(isHydrostatic=self.wvt.isHydrostatic), "AbsTol",1e-7,"RelTol",1e-7);
            self.verifyEqual(self.wvt.totalEnstrophy,soln.depthIntegratedTotalEnstrophy, "AbsTol",1e-7,"RelTol",1e-7);
        end

        function testAddWaveModes(self)
            self.wvt.removeAll();

            kMode = 3; lMode = 5; j=3; phi = pi*0.3; u=0.05; sign = +1;
            self.wvt.setWaveModes(kMode=kMode,lMode=lMode,j=j,phi=phi,u=u,sign=sign);
            soln1 = self.wvt.waveComponent.internalGravityWaveSolution(kMode,lMode,j,u,phi,sign,amplitudeIsMaxU=1,t=self.wvt.t);

            kMode = 4; lMode = -1; j=2; phi = pi*0.1; u=0.05; sign = +1;
            self.wvt.addWaveModes(kMode=kMode,lMode=lMode,j=j,phi=phi,u=u,sign=sign);
            soln2 = self.wvt.waveComponent.internalGravityWaveSolution(kMode,lMode,j,u,phi,sign,amplitudeIsMaxU=1,t=self.wvt.t);

            self.wvt.t = 86400;
            args = {self.wvt.X,self.wvt.Y,self.wvt.Z,self.wvt.t};
            self.verifyThat(self.wvt.u,IsSameSolutionAs(soln1.u(args{:})+soln2.u(args{:})),'u');
            self.verifyThat(self.wvt.v,IsSameSolutionAs(soln1.v(args{:})+soln2.v(args{:})),'v');
            self.verifyThat(self.wvt.w,IsSameSolutionAs(soln1.w(args{:})+soln2.w(args{:})),'w');
            self.verifyThat(self.wvt.eta,IsSameSolutionAs(soln1.eta(args{:})+soln2.eta(args{:})),'eta');
            self.verifyThat(self.wvt.p,IsSameSolutionAs(soln1.p(args{:})+soln2.p(args{:})),'p');
            self.verifyThat(self.wvt.qgpv,IsSameSolutionAs(soln1.qgpv(args{:})+soln2.qgpv(args{:})),'qgpv');

            self.verifyEqual(self.wvt.totalEnergy,soln1.depthIntegratedTotalEnergy(isHydrostatic=self.wvt.isHydrostatic)+soln2.depthIntegratedTotalEnergy(isHydrostatic=self.wvt.isHydrostatic), "AbsTol",1e-7,"RelTol",1e-3);
            self.verifyEqual(self.wvt.totalEnstrophy,soln1.depthIntegratedTotalEnstrophy+soln2.depthIntegratedTotalEnstrophy, "AbsTol",1e-7,"RelTol",1e-7);
        end

    end

end