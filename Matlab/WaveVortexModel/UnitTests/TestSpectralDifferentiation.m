classdef TestSpectralDifferentiation < matlab.unittest.TestCase
    properties
        wvt
    end

    properties (ClassSetupParameter)
        % transform = {'constant','hydrostatic','boussinesq'};
        transform = {'constant'};
    end

    methods (TestClassSetup)
        function classSetup(testCase,transform)
            switch transform
                case 'constant'
                    testCase.wvt = WVTransformConstantStratification([1, 10, 4], [16, 8, 17]);
                case 'hydrostatic'
                    testCase.wvt = WVTransformHydrostatic([1, 10, 4], [16, 8, 17], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
                case 'boussinesq'
                    testCase.wvt = WVTransformBoussinesq([1, 10, 4], [16, 8, 17], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));
            end
        end
    end

    properties (TestParameter)
        derivative = struct('first',1,'second',2,'third',3,'fourth',4);
    end

    methods (Test)
        function testDiffX(testCase,derivative)
            [X,Y,Z] = testCase.wvt.xyzGrid;
            Lx = testCase.wvt.Lx;
            Ly = testCase.wvt.Ly;

            n = 1;
            for m = (1:floor(testCase.wvt.Nx/2))
                kx = 2*pi*m/Lx;
                ky = 2*pi*n/Ly;
                f = cos(kx*X) .* cos(ky*Y);
                switch derivative
                    case 1
                        Df_analytical = -kx*sin(kx*X).*cos(ky*Y);
                    case 2
                        Df_analytical = -(kx^2)*cos(kx*X).*cos(ky*Y);
                    case 3
                        Df_analytical = (kx^3)*sin(kx*X).*cos(ky*Y);
                    case 4
                        Df_analytical = (kx^4)*cos(kx*X).*cos(ky*Y);
                end

                testCase.verifyEqual(testCase.wvt.diffX(f,derivative),Df_analytical, ...
                    "AbsTol",1e-7,"RelTol",1e-7);
            end
        end
    end

end