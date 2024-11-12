classdef TestNonlinearFlux < matlab.unittest.TestCase
    properties
        wvt_
    end

    properties (ClassSetupParameter)
        Lxyz = struct('Lxyz',[4e3, 4e3, 2e3]);
        % Nxyz = struct('Nx8Ny8Nz5',[8 8 5]);
        Nxyz = struct('Nx16Ny16Nz9',[16 16 9]);
        %transform = {'constant','hydrostatic','boussinesq'};
        transform = {'constant-hydrostatic','constant-boussinesq','hydrostatic','boussinesq'};
    end

    methods (TestClassSetup)
        function classSetup(testCase,Lxyz,Nxyz,transform)
            switch transform
                case 'constant-hydrostatic'
                    testCase.wvt_ = WVTransformConstantStratification(Lxyz, Nxyz, isHydrostatic=1, shouldAntialias=0);
                case 'constant-boussinesq'
                    testCase.wvt_ = WVTransformConstantStratification(Lxyz, Nxyz, isHydrostatic=0, shouldAntialias=0);
                    testCase.wvt_.nonlinearFluxOperation = WVNonlinearFluxNonhydrostatic(testCase.wvt_);
                case 'hydrostatic'
                    testCase.wvt_ = WVTransformHydrostatic(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)),shouldAntialias=0);
                case 'boussinesq'
                    testCase.wvt_ = WVTransformBoussinesq(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)),shouldAntialias=0);
            end
        end
    end

    methods (Test)
        function testNonlinearFlux(self)
            wvt = self.wvt_;
            wvt.initWithRandomFlow();
            spatialFlux = WVNonlinearFluxSpatial(wvt);
            standardFlux = WVNonlinearFlux(wvt);

            wvt.t = 6000;
            [SpatialFp,SpatialFm,SpatialF0] = spatialFlux.compute(wvt);
            [StandardFp,StandardFm,StandardF0] = standardFlux.compute(wvt);

            self.verifyEqual(StandardFp,SpatialFp, "AbsTol",1e-7,"RelTol",1e-7);
            self.verifyEqual(StandardFm,SpatialFm, "AbsTol",1e-7,"RelTol",1e-7);
            self.verifyEqual(StandardF0,SpatialF0, "AbsTol",1e-7,"RelTol",1e-7);
        end

        function testEnergyFluxConservation(self)
            wvt = self.wvt_;
            wvt.initWithRandomFlow('wave',uvMax=0.1);
            
            % We are careful to *not* initialize in an anti-aliased
            % configuration, and then only energize modes that will not
            % alias. This ensures that energy can be conserved.
            antialiasMask = zeros(wvt.spectralMatrixSize);
            antialiasMask(wvt.Kh > 2*max(abs(wvt.k))/3) = 1;
            antialiasMask(wvt.J > 2*max(abs(wvt.j))/3) = 1;
            antialiasMask = logical(antialiasMask);

            wvt.Ap(antialiasMask) = 0;
            wvt.Am(antialiasMask) = 0;
            wvt.A0(antialiasMask) = 0;

            [Fp,Fm,F0] = wvt.nonlinearFlux();
            [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0,deltaT=0);
            totalEnergyFlux = sum(Ep(:))+sum(Em(:))+sum(E0(:));
            self.verifyEqual(totalEnergyFlux,0, "AbsTol",1e-15,"RelTol",1e-7);
        end

        % function testNonlinearWaveTriad(self)
        %     wvt = self.wvt_;
        % end

    end

end