classdef TestSpectralDifferentiationXY < matlab.unittest.TestCase
    properties
        wvt
    end

    properties (ClassSetupParameter)
        %transform = {'constant','hydrostatic','boussinesq'};
        % Lxyz = struct('Lxyz',[1 10 4]);
        Nxyz = struct('Nx16Ny16Nz9',[16 16 9]);
        % transform = {'hydrostatic'};
        Lxyz = struct('Lxyz',[1000, 500, 500]);
        % Nxyz = struct('Nx32N16Nz17',[32 16 17]);
        transform = {'hydrostatic'};
    end

    methods (TestClassSetup)
        function classSetup(testCase,Lxyz,Nxyz,transform)
            switch transform
                case 'constant'
                    testCase.wvt = WVTransformConstantStratification(Lxyz, Nxyz,shouldAntialias=false);
                case 'hydrostatic'
                    testCase.wvt = WVTransformHydrostatic(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)),shouldAntialias=false);
                case 'boussinesq'
                    testCase.wvt = WVTransformBoussinesq(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)),shouldAntialias=false);
            end
        end
    end

    methods (TestParameterDefinition,Static)
        function [k_n,l_n] = initializeProperty(Lxyz,Nxyz,transform)
            % If you want to dynamically adjust the test parameters, you
            % have to do it here.
            for i=1:(floor(Nxyz(1)/2)-1) % Note that we are specifically avoiding testing the Nyquist which is not fully resolved.
                k_n.(sprintf('k_%d',i)) = i;
            end
            for i=1:(floor(Nxyz(2)/2)-1) % Note that we are specifically avoiding testing the Nyquist which is not fully resolved.
                l_n.(sprintf('l_%d',i)) = i;
            end
        end
    end

    properties (TestParameter)
        derivative = struct('first',1,'second',2,'third',3,'fourth',4);
        % derivative = struct('second',2,'third',3,'fourth',4);
        k_n
        l_n
    end

    methods (Test)
        function testDiffX(testCase,derivative,k_n,l_n)
            [X,Y,Z] = testCase.wvt.xyzGrid;
            Lx = testCase.wvt.Lx;
            Ly = testCase.wvt.Ly;
            phix = 2*pi*rand(1);
            phiy = 2*pi*rand(1);

            kx = 2*pi*k_n/Lx;
            ky = 2*pi*l_n/Ly;
            f = cos(kx*X+phix) .* cos(ky*Y+phiy);
            switch derivative
                case 1
                    Df_analytical = -kx*sin(kx*X+phix).*cos(ky*Y+phiy);
                case 2
                    Df_analytical = -(kx^2)*cos(kx*X+phix).*cos(ky*Y+phiy);
                case 3
                    Df_analytical = (kx^3)*sin(kx*X+phix).*cos(ky*Y+phiy);
                case 4
                    Df_analytical = (kx^4)*cos(kx*X+phix).*cos(ky*Y+phiy);
            end

            testCase.verifyEqual(testCase.wvt.diffX(f,n=derivative),Df_analytical, "AbsTol",1e-7,"RelTol",1e-7);
        end

        function testDiffY(testCase,derivative,k_n,l_n)
            [X,Y,Z] = testCase.wvt.xyzGrid;
            Lx = testCase.wvt.Lx;
            Ly = testCase.wvt.Ly;
            phix = 2*pi*rand(1);
            phiy = 2*pi*rand(1);


            kx = 2*pi*k_n/Lx;
            ky = 2*pi*l_n/Ly;
            f = cos(kx*X+phix) .* cos(ky*Y+phiy);
            switch derivative
                case 1
                    Df_analytical = -ky*sin(ky*Y+phiy).*cos(kx*X+phix);
                case 2
                    Df_analytical = -(ky^2)*cos(ky*Y+phiy).*cos(kx*X+phix);
                case 3
                    Df_analytical = (ky^3)*sin(ky*Y+phiy).*cos(kx*X+phix);
                case 4
                    Df_analytical = (ky^4)*cos(ky*Y+phiy).*cos(kx*X+phix);
            end

            testCase.verifyEqual(testCase.wvt.diffY(f,n=derivative),Df_analytical, "AbsTol",1e-7,"RelTol",1e-7);

        end
    end

end