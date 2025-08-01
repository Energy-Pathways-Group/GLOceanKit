classdef TestGeostrophicMethods < matlab.unittest.TestCase
    properties
        wvt
        solutionGroup
    end

    properties (ClassSetupParameter)
        Lxyz = struct('Lxyz',[750e3, 750e3, 1300]);
        Nxyz = struct('Nx64Ny64Nz30',[64 64 40]);
        % Nxyz = struct('Nx16Ny16Nz5',[16 16 5]);
        % transform = {'constant','hydrostatic','boussinesq'};
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
            % testCase.wvt.addOperation(testCase.wvt.operationForDynamicalVariable('u','v','eta','w',flowComponent=testCase.wvt.flowComponent('geostrophic')));
            testCase.solutionGroup = WVGeostrophicComponent(testCase.wvt);
        end
    end

    methods (Test)
        function testRemoveAllGeostrophicMotions(self)
            % In this test we intialize with a random flow state, confirm
            % that both total energy and geostrophic energy are present,
            % remove all the geostrophic energy, and then confirm that there
            % is no geostrophic energy remaining, and that the total energy is
            % the same as the initial total, minus the initial geostrophic.
            self.wvt.removeAll();

            % currently we need to make this small-ish to avoid density
            % overturns.
            self.wvt.initWithRandomFlow(uvMax=0.025);

            initialTotalEnergy = self.wvt.totalEnergy;
            initialGeostrophicEnergy = self.wvt.geostrophicEnergy;
            self.verifyGreaterThan(initialTotalEnergy,0.0);
            self.verifyGreaterThan(initialGeostrophicEnergy,0.0);

            self.wvt.removeAllGeostrophicMotions();
            finalTotalEnergy = self.wvt.totalEnergy;
            finalGeostrophicEnergy = self.wvt.geostrophicEnergy;

            self.verifyEqual(finalGeostrophicEnergy,0.0);
            self.verifyEqual(finalTotalEnergy,initialTotalEnergy-initialGeostrophicEnergy, "AbsTol",1e-7,"RelTol",1e-7);
        end

        function testMeanPressureViolation(self)
            % zero out our pressure correction, and confirm we throw an
            % error
            self.wvt.removeAll();

            Le = 120e3;
            He = 200;
            x0 = (1/2)*max(self.wvt.x); y0=max(self.wvt.y)/2;

            U = 0.10; % m/s
            psibar = @(z) (pi*Le*Le/(self.wvt.Lx*self.wvt.Ly))*U*(Le/sqrt(2))*exp(1/2)*exp(-(z/He/sqrt(2)).^2 );
            psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2 -(z/He/sqrt(2)).^2 ) - 0*psibar(z);
            self.verifyError(@() self.wvt.initWithGeostrophicStreamfunction(psi),'WVTransform:MeanPressureViolation' );
        end

        function testMeanPressureNonViolation(self)
            % confirm no error with the pressure correction
            self.wvt.removeAll();

            Le = 120e3;
            He = 300;
            x0 = (1/2)*max(self.wvt.x); y0=max(self.wvt.y)/2;

            U = 0.10; % m/s
            psibar = @(z) (pi*Le*Le/(self.wvt.Lx*self.wvt.Ly))*U*(Le/sqrt(2))*exp(1/2)*exp(-(z/He/sqrt(2)).^2 );
            psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2 -(z/He/sqrt(2)).^2 ) - psibar(z);
            self.verifyWarningFree(@() self.wvt.initWithGeostrophicStreamfunction(psi));
        end

        function testDensityBoundsViolation(self)
            % confirm an error when the amplitude is too large
            self.wvt.removeAll();

            Le = 120e3;
            He = 200;
            x0 = (1/2)*max(self.wvt.x); y0=max(self.wvt.y)/2;

            U = 0.20; % m/s
            psibar = @(z) (pi*Le*Le/(self.wvt.Lx*self.wvt.Ly))*U*(Le/sqrt(2))*exp(1/2)*exp(-(z/He/sqrt(2)).^2 );
            psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2 -(z/He/sqrt(2)).^2 ) - psibar(z);
            self.verifyError(@() self.wvt.initWithGeostrophicStreamfunction(psi),'WVTransform:DensityBoundsViolation' );

            self.wvt.removeAll();
            U = -1.00; % m/s
            psibar = @(z) (pi*Le*Le/(self.wvt.Lx*self.wvt.Ly))*U*(Le/sqrt(2))*exp(1/2)*exp(-(z/He/sqrt(2)).^2 );
            psi = @(x,y,z) U*(Le/sqrt(2))*exp(1/2)*exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2 -(z/He/sqrt(2)).^2 ) - psibar(z);
            self.verifyError(@() self.wvt.initWithGeostrophicStreamfunction(psi),'WVTransform:DensityBoundsViolation' );
        end

        function testInitWithGeostrophicStreamfunction(self)
            % In this test we expect all existing motions to be removed
            % when we initialize.
            %
            % VERY IMPORTANT: We do *not* get back exactly what we put in
            % because we are de-aliasing the signal. Hence, this experiment
            % was designed assuming Nz=30 (and thus Nj=20).
            self.wvt.removeAll();

            Le = 120e3;
            He = 300;
            x0 = (1/2)*max(self.wvt.x); y0=max(self.wvt.y)/2;
            rho0 = self.wvt.rho0;
            f = self.wvt.f;
            g = self.wvt.g;

            U = 0.10; % m/s
            A = U*(Le/sqrt(2))*exp(1/2);
            H = @(z) exp(-(z/He/sqrt(2)).^2 );
            F = @(x,y) exp(-((x-x0)/Le).^2 -((y-y0)/Le).^2);
            psibar = (pi*Le*Le/(self.wvt.Lx*self.wvt.Ly));
            psi = @(x,y,z) A*H(z).*(F(x,y) - psibar);
            u = @(x,y,z) A*(2/Le/Le)*(y-y0).*H(z).*F(x,y);
            v = @(x,y,z) -A*(2/Le/Le)*(x-x0).*H(z).*F(x,y);
            rho_e = @(x,y,z) (rho0*f/g)*(1/He/He)*A*z.*H(z).*(F(x,y) - psibar);
            % Populate the flow field with junk...
            % self.wvt.initWithRandomFlow();

            % call our initWithInertialMotions method
            self.wvt.initWithGeostrophicStreamfunction(psi);

            % now verify that only inertial oscillations are part of the
            % solution.
            self.verifyThat(u(self.wvt.X,self.wvt.Y,self.wvt.Z),IsSameSolutionAs(self.wvt.u,relTol=1e-3),'u');
            self.verifyThat(v(self.wvt.X,self.wvt.Y,self.wvt.Z),IsSameSolutionAs(self.wvt.v,relTol=1e-3),'v');
            self.verifyThat(rho_e(self.wvt.X,self.wvt.Y,self.wvt.Z),IsSameSolutionAs(self.wvt.rho_e,relTol=1e-3),'rho_e');
            self.verifyThat(psi(self.wvt.X,self.wvt.Y,self.wvt.Z),IsSameSolutionAs(self.wvt.psi,relTol=1e-3),'psi');
        end

        function testSetGeostrophicModes(self)
            self.wvt.removeAll();

            kMode = 3; lMode = 5; j=3; phi = pi*0.3; u=0.1;

            self.wvt.setGeostrophicModes(kMode=kMode,lMode=lMode,j=j,phi=phi,u=u);

            soln = self.wvt.geostrophicComponent.geostrophicSolution(kMode,lMode,j,u,phi,amplitudeIsMaxU=1);

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

        function testAddGeostrophicModes(self)
            self.wvt.removeAll();

            kMode = 3; lMode = 5; j=3; phi = pi*0.3; u=0.1;
            self.wvt.setGeostrophicModes(kMode=kMode,lMode=lMode,j=j,phi=phi,u=u);
            soln1 = self.wvt.geostrophicComponent.geostrophicSolution(kMode,lMode,j,u,phi,amplitudeIsMaxU=1);

            kMode = 4; lMode = 1; j=2; phi = pi*0.1; u=0.05;
            self.wvt.addGeostrophicModes(kMode=kMode,lMode=lMode,j=j,phi=phi,u=u);
            soln2 = self.wvt.geostrophicComponent.geostrophicSolution(kMode,lMode,j,u,phi,amplitudeIsMaxU=1);

            self.wvt.t = 86400;
            args = {self.wvt.X,self.wvt.Y,self.wvt.Z,self.wvt.t};
            self.verifyThat(self.wvt.u,IsSameSolutionAs(soln1.u(args{:})+soln2.u(args{:})),'u');
            self.verifyThat(self.wvt.v,IsSameSolutionAs(soln1.v(args{:})+soln2.v(args{:})),'v');
            self.verifyThat(self.wvt.w,IsSameSolutionAs(soln1.w(args{:})+soln2.w(args{:})),'w');
            self.verifyThat(self.wvt.eta,IsSameSolutionAs(soln1.eta(args{:})+soln2.eta(args{:})),'eta');
            self.verifyThat(self.wvt.p,IsSameSolutionAs(soln1.p(args{:})+soln2.p(args{:})),'p');
            self.verifyThat(self.wvt.qgpv,IsSameSolutionAs(soln1.qgpv(args{:})+soln2.qgpv(args{:})),'qgpv');

            self.verifyEqual(self.wvt.totalEnergy,soln1.depthIntegratedTotalEnergy(isHydrostatic=self.wvt.isHydrostatic)+soln2.depthIntegratedTotalEnergy(isHydrostatic=self.wvt.isHydrostatic), "AbsTol",1e-7,"RelTol",1e-3);
            self.verifyEqual(self.wvt.totalEnstrophy,soln1.depthIntegratedTotalEnstrophy+soln2.depthIntegratedTotalEnstrophy, "AbsTol",1e-7,"RelTol",1e-3);
        end

    end

end