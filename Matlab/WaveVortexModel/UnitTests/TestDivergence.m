classdef TestDivergence < matlab.unittest.TestCase

    methods (Test)
        function testDivergence(self)
            Lxyz = [1000, 500, 500];
            Nxyz = [16 8 9];
            latitude = 33;
            % wvt = WVTransformConstantStratification(Lxyz, Nxyz, latitude=latitude, isHydrostatic=0,shouldAntialias=0);
            wvt = WVTransformHydrostatic(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)),latitude=latitude);

            wvt.initWithRandomFlow(uvMax=0.1);

            delta = -wvt.diffX(wvt.u) - wvt.diffY(wvt.v);
            delta_bar = wvt.transformFromSpatialDomainWithFourier(delta);
            delta_hat = wvt.transformFromSpatialDomainWithFg(delta_bar);
            w = wvt.transformToSpatialDomainWithG(A0=wvt.h_0 .* delta_hat);
            
            self.verifyEqual(w,wvt.w, "AbsTol",1e-7,"RelTol",1e-7);
        end

        function testForwardAndInverseWVTrasform(self)
            Lxyz = [1000, 500, 500];
            Nxyz = [16 8 9];
            latitude = 33;
            % wvt = WVTransformConstantStratification(Lxyz, Nxyz, latitude=latitude, isHydrostatic=0,shouldAntialias=0);
            wvt = WVTransformHydrostatic(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)),latitude=latitude);
            
            wvt.initWithRandomFlow(uvMax=0.1);

            [ApNL,AmNL,A0NL] = wvt.transformUVEtaToWaveVortex(wvt.u,wvt.v,wvt.eta);
            [u,v,w,eta] = wvt.transformWaveVortexToUVWEta(ApNL,AmNL,A0NL);
            self.verifyEqual(u,wvt.u, "AbsTol",1e-7,"RelTol",1e-7);
            self.verifyEqual(v,wvt.v, "AbsTol",1e-7,"RelTol",1e-7);
            self.verifyEqual(w,wvt.w, "AbsTol",1e-7,"RelTol",1e-7);
            self.verifyEqual(eta,wvt.eta, "AbsTol",1e-7,"RelTol",1e-7);
        end

    end

end