classdef TestSpectralDifferentiationZ < matlab.unittest.TestCase
    properties
        wvt
    end

    properties (ClassSetupParameter)
        % transform = {'constant','hydrostatic','boussinesq'};
        % Lxyz = struct('Lxyz',[1 10 4]);
        % Nxyz = struct('Nx16Ny16Nz9',[16 16 9]);
        % transform = {'constant','hydrostatic','boussinesq'};
        Lxyz = struct('Lxyz',[1000, 500, 500]);
        Nxyz = struct('Nx32N16Nz17',[32 16 17]);
        transform = {'hydrostatic'};
    end

    methods (TestClassSetup)
        function classSetup(testCase,Lxyz,Nxyz,transform)
            switch transform
                case 'constant'
                    testCase.wvt = WVTransformConstantStratification(Lxyz, Nxyz, shouldAntialias=0);
                case 'hydrostatic'
                    testCase.wvt = WVTransformHydrostatic(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)),shouldAntialias=0);
                case 'boussinesq'
                    testCase.wvt = WVTransformBoussinesq(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)),shouldAntialias=0);
            end
        end
    end

    methods (TestParameterDefinition,Static)
        function [k_n,l_n,m_n] = initializeProperty(Lxyz,Nxyz,transform)
            % If you want to dynamically adjust the test parameters, you
            % have to do it here.
            for i=1:floor(Nxyz(1)/2)
                k_n.(sprintf('k_%d',i)) = i;
            end
            for i=1:floor(Nxyz(2)/2)
                l_n.(sprintf('l_%d',i)) = i;
            end
            for i=0:(Nxyz(3)-1)
                m_n.(sprintf('m_%d',i)) = i;
            end
        end
    end

    properties (TestParameter)
        derivative = struct('first',1,'second',2,'third',3,'fourth',4);
        %derivative = struct('first',1);
        k_n
        l_n
        m_n
    end

    methods (Test)
        function testDiffZF(self,derivative,k_n,m_n)
            [X,Y,Z] = self.wvt.xyzGrid;

            Lx = self.wvt.Lx;
            Ly = self.wvt.Ly;
            df = 1/((self.wvt.Nz-1)*(self.wvt.z(2)-self.wvt.z(1)));

            phix = 2*pi*rand(1);
            phiy = 2*pi*rand(1);
            n = 1;

            kx = 2*pi*k_n/Lx;
            ky = 2*pi*n/Ly;
            kz = pi*df*m_n;
            f = cos(kx*X+phix) .* cos(ky*Y+phiy) .* cos(kz*Z);
            switch derivative
                case 1
                    Df_analytical = -kz*sin(kz*Z).*cos(kx*X+phix).*cos(ky*Y+phiy);
                case 2
                    Df_analytical = -(kz^2)*cos(kz*Z).*cos(kx*X+phix).*cos(ky*Y+phiy);
                case 3
                    Df_analytical = (kz^3)*sin(kz*Z).*cos(kx*X+phix).*cos(ky*Y+phiy);
                case 4
                    Df_analytical = (kz^4)*cos(kz*Z).*cos(kx*X+phix).*cos(ky*Y+phiy);
            end

            self.verifyEqual(self.wvt.diffZF(f,n=derivative),Df_analytical, "AbsTol",1e-7,"RelTol",1e-7);

        end

        function testDiffZG(self,derivative,k_n,m_n)
            [X,Y,Z] = self.wvt.xyzGrid;

            Lx = self.wvt.Lx;
            Ly = self.wvt.Ly;
            df = 1/((self.wvt.Nz-1)*(self.wvt.z(2)-self.wvt.z(1)));

            phix = 2*pi*rand(1);
            phiy = 2*pi*rand(1);
            n = 1;

            kx = 2*pi*k_n/Lx;
            ky = 2*pi*n/Ly;
            kz = pi*df*m_n;
            f = cos(kx*X+phix) .* cos(ky*Y+phiy) .* sin(kz*Z);
            switch derivative
                case 1
                    Df_analytical = kz*cos(kz*Z).*cos(kx*X+phix).*cos(ky*Y+phiy);
                case 2
                    Df_analytical = -(kz^2)*sin(kz*Z).*cos(kx*X+phix).*cos(ky*Y+phiy);
                case 3
                    Df_analytical = -(kz^3)*cos(kz*Z).*cos(kx*X+phix).*cos(ky*Y+phiy);
                case 4
                    Df_analytical = (kz^4)*sin(kz*Z).*cos(kx*X+phix).*cos(ky*Y+phiy);
            end

            if m_n==16 && (derivative == 1 || derivative == 3)
                % Nquist, cosine of m=8 is not resolved.
                self.verifyNotEqual(self.wvt.diffZG(f,n=derivative),Df_analytical);
            else
                self.verifyEqual(self.wvt.diffZG(f,n=derivative),Df_analytical, "AbsTol",1e-7,"RelTol",1e-7);
            end

        end


    end

end