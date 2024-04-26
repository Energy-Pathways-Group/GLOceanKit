classdef TestNonlinearFlux < matlab.unittest.TestCase
    properties
        wvt
    end

    properties (ClassSetupParameter)
        Lxyz = struct('Lxyz',[15e3, 15e3, 4e3]);
        % Nxyz = struct('Nx8Ny8Nz5',[8 8 5]);
        Nxyz = struct('Nx16Ny16Nz9',[16 16 9]);
        %transform = {'constant','hydrostatic','boussinesq'};
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
        end
    end

    methods (Test)
        function testNonlinearFlux(self)
            self.wvt.initWithRandomFlow();
            spatialFlux = WVNonlinearFluxSpatial(self.wvt);
            standardFlux = WVNonlinearFlux(self.wvt,shouldAntialias=0);

            [SpatialFp,SpatialFm,SpatialF0] = spatialFlux.compute(self.wvt);
            [StandardFp,StandardFm,StandardF0] = standardFlux.compute(self.wvt);

            self.verifyEqual(StandardFp,SpatialFp, "AbsTol",1e-7,"RelTol",1e-7);
            self.verifyEqual(StandardFm,SpatialFm, "AbsTol",1e-7,"RelTol",1e-7);
            self.verifyEqual(StandardF0,SpatialF0, "AbsTol",1e-7,"RelTol",1e-7);
        end

    end

end