classdef TestFourierTransformXY < matlab.unittest.TestCase
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
                    testCase.wvt = WVTransformConstantStratification(Lxyz, Nxyz, shouldAntialias=false);
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
            for i=0:(floor(Nxyz(1)/2)-1) % Note that we are specifically avoiding testing the Nyquist which is not fully resolved.
                k_n.(sprintf('k_%d',i)) = i;
            end
            for i=0:(floor(Nxyz(2)/2)-1) % Note that we are specifically avoiding testing the Nyquist which is not fully resolved.
                l_n.(sprintf('l_%d',i)) = i;
            end
        end
    end

    properties (TestParameter)
        k_n
        l_n
    end

    methods (Test)
        function testForwardBackwardTransform(testCase,k_n,l_n)
            [X,Y,Z] = testCase.wvt.xyzGrid;
            Lx = testCase.wvt.Lx;
            Ly = testCase.wvt.Ly;
            phix = 2*pi*rand(1);
            phiy = 2*pi*rand(1);

            kx = 2*pi*k_n/Lx;
            ky = 2*pi*l_n/Ly;
            f = cos(kx*X+phix) .* cos(ky*Y+phiy);

            f_bar = testCase.wvt.transformFromSpatialDomainWithFourier(f);
            f_back = testCase.wvt.transformToSpatialDomainWithFourier(f_bar);

            testCase.verifyEqual(f,f_back, "AbsTol",1e-7,"RelTol",1e-7);
        end
    end

end