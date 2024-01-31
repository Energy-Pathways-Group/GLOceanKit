classdef TestSpectralDifferentiationZ < matlab.unittest.TestCase
    properties
        wvt
    end

    properties (ClassSetupParameter)
        % transform = {'constant','hydrostatic','boussinesq'};
        transform = {'boussinesq'};
    end

    methods (TestClassSetup)
        function classSetup(testCase,transform)
            switch transform
                case 'constant'
                    testCase.wvt = WVTransformConstantStratification([1, 10, 4], [16 16 9]);
                case 'hydrostatic'
                    testCase.wvt = WVTransformHydrostatic([1, 10, 4], [16 16 9], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
                case 'boussinesq'
                    testCase.wvt = WVTransformBoussinesq([1, 10, 4], [16 16 9], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
            end
        end
    end

    properties (TestParameter)
        % derivative = struct('first',1,'second',2,'third',3,'fourth',4);
        derivative = struct('first',1);
        k_n = struct('k_1',1,'k_2',2,'k_3',3,'k_4',4,'k_5',5,'k_6',6,'k_7',7,'k_8',8)
        l_n = struct('l_1',1,'l_2',2,'l_3',3,'l_4',4,'l_5',5,'l_6',6,'l_7',7,'l_8',8)
        m_n = struct('m_0',0,'m_1',1,'m_2',2,'m_3',3,'m_4',4,'m_5',5,'m_6',6,'m_7',7,'m_8',8)
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

            self.verifyEqual(self.wvt.diffZF(f,derivative),Df_analytical, "AbsTol",1e-7,"RelTol",1e-7);

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

            if m_n==8 && (derivative == 1 || derivative == 3)
                % Nquist, cosine of m=8 is not resolved.
                self.verifyNotEqual(self.wvt.diffZG(f,derivative),Df_analytical);
            else
                self.verifyEqual(self.wvt.diffZG(f,derivative),Df_analytical, "AbsTol",1e-7,"RelTol",1e-7);
            end

        end


    end

end