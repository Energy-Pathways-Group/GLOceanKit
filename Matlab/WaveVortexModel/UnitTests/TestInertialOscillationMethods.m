classdef TestInertialOscillationMethods < matlab.unittest.TestCase
    properties
        wvt
        solutionGroup
    end

    properties (ClassSetupParameter)
        Lxyz = struct('Lxyz',[15e3, 15e3, 1300]);
        Nxyz = struct('Nx8Ny8Nz5',[8 8 10]);
        % Nxyz = struct('Nx16Ny16Nz5',[16 16 5]);
        % transform = {'constant','hydrostatic','boussinesq'};
        transform = {'constant'};
    end

    methods (TestClassSetup)
        function classSetup(testCase,Lxyz,Nxyz,transform)
            switch transform
                case 'constant'
                    testCase.wvt = WVTransformConstantStratification(Lxyz, Nxyz);
                case 'hydrostatic'
                    testCase.wvt = WVTransformHydrostatic(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)),shouldAntialias=0);
                case 'boussinesq'
                    testCase.wvt = WVTransformBoussinesq(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
            end
            testCase.solutionGroup = WVInertialOscillationComponent(testCase.wvt);
        end
    end

    methods (Test)
        function testRemoveAllInertialMotions(self)
            % In this test we intialize with a random flow state, confirm
            % that both total energy and inertial energy are present,
            % remove all the inertial energy, and then confirm that there
            % is no inertial energy remaining, and that the total energy is
            % the same as the initial total, minus the initial inertial.
            self.wvt.initWithRandomFlow();

            initialTotalEnergy = self.wvt.totalEnergy;
            initialInertialEnergy = self.wvt.inertialEnergy;
            self.verifyGreaterThan(initialTotalEnergy,0.0);
            self.verifyGreaterThan(initialInertialEnergy,0.0);

            self.wvt.removeAllInertialMotions();
            finalTotalEnergy = self.wvt.totalEnergy;
            finalInertialEnergy = self.wvt.inertialEnergy;

            self.verifyEqual(finalInertialEnergy,0.0);
            self.verifyEqual(finalTotalEnergy,initialTotalEnergy-initialInertialEnergy, "AbsTol",1e-7,"RelTol",1e-7);
        end

        function testInitWithInertialMotions(self)
            U_io = 0.2;
            Ld = self.wvt.Lz/5;
            theta = 0;
            u_NIO = @(z) U_io*cos(theta)*exp((z/Ld));
            v_NIO = @(z) U_io*sin(theta)*exp((z/Ld));

            self.wvt.addOperation(self.wvt.operationForDynamicalVariable('u','v','eta','w',flowComponent=self.wvt.flowComponent('inertial')));

            self.wvt.initWithRandomFlow();
            % self.wvt.initWithInertialMotions(u_NIO,v_NIO);
            self.wvt.initWithUVEta(u_NIO(self.wvt.Z),v_NIO(self.wvt.Z),zeros(self.wvt.spatialMatrixSize));

            self.verifyThat(u_NIO(self.wvt.Z),IsSameSolutionAs(self.wvt.u),'u_tot');
            self.verifyThat(v_NIO(self.wvt.Z),IsSameSolutionAs(self.wvt.v),'v_tot');
        end

        % function testSetInertialMotions(self)
        %     U_io = 0.2;
        %     Ld = self.wvt.Lz/5;
        %     theta = 0;
        %     u_NIO = @(z) U_io*cos(theta)*exp((z/Ld));
        %     v_NIO = @(z) U_io*sin(theta)*exp((z/Ld));
        % 
        %     self.wvt.addOperation(self.wvt.operationForDynamicalVariable('u','v','eta','w',flowComponent=self.wvt.flowComponent('inertial')));
        % 
        %     self.wvt.initWithRandomFlow();
        % 
        %     initialTotalEnergy = self.wvt.totalEnergy;
        %     initialInertialEnergy = self.wvt.inertialEnergy;
        % 
        %     self.wvt.setInertialMotions(u_NIO,v_NIO);
        % 
        %     finalTotalEnergy = self.wvt.totalEnergy;
        %     finalInertialEnergy = self.wvt.inertialEnergy;
        % 
        %     self.verifyThat(u_NIO(self.wvt.Z),IsSameSolutionAs(self.wvt.u_io),'u_io');
        %     self.verifyThat(v_NIO(self.wvt.Z),IsSameSolutionAs(self.wvt.v_io),'v_io');
        %     self.verifyEqual(finalTotalEnergy-finalInertialEnergy,initialTotalEnergy-initialInertialEnergy, "AbsTol",1e-7,"RelTol",1e-7);
        % end

    end

end