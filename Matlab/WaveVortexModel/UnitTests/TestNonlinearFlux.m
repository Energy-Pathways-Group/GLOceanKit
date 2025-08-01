classdef TestNonlinearFlux < matlab.unittest.TestCase
    properties
        wvt_
    end

    properties (ClassSetupParameter)
        Lxyz = struct('Lxyz',[4e3, 4e3, 2e3]);
        % Nxyz = struct('Nx8Ny8Nz5',[8 8 5]);
        Nxyz = struct('Nx16Ny16Nz9',[16 16 9]);
        %transform = {'constant','hydrostatic','boussinesq'};
        % transform = {'constant-hydrostatic','constant-boussinesq','hydrostatic','boussinesq'};
        % transform = {'constant-hydrostatic','constant-boussinesq'};
        % transform = {'boussinesq'};
        % transform = {'hydrostatic'};
        transform = {'hydrostatic-exp'};
    end

    methods (TestClassSetup)
        function classSetup(testCase,Lxyz,Nxyz,transform)
            switch transform
                case 'constant-hydrostatic'
                    testCase.wvt_ = WVTransformConstantStratification(Lxyz, Nxyz, isHydrostatic=true, shouldAntialias=0);
                case 'constant-boussinesq'
                    testCase.wvt_ = WVTransformConstantStratification(Lxyz, Nxyz,shouldAntialias=0);
                case 'hydrostatic'
                    testCase.wvt_ = WVTransformHydrostatic(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)),shouldAntialias=false);
                case 'hydrostatic-exp'
                    N0 = 3*2*pi/3600; % buoyancy frequency at the surface, radians/seconds
                    L_gm = 1300; % thermocline exponential scale, meters
                    N2 = @(z) N0*N0*exp(2*z/L_gm);
                    testCase.wvt_ = WVTransformHydrostatic(Lxyz, Nxyz, N2=N2,shouldAntialias=false);
                case 'boussinesq'
                    testCase.wvt_ = WVTransformBoussinesq(Lxyz, Nxyz, N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)),shouldAntialias=0);
            end
        end
    end

    methods (Test)
        % function testNonlinearFlux(self)
        %     wvt = self.wvt_;
        %     wvt.initWithRandomFlow();
        %     spatialFlux = WVNonlinearFluxSpatial(wvt);
        %     standardFlux = WVNonlinearFlux(wvt);
        % 
        %     wvt.t = 6000;
        %     [SpatialFp,SpatialFm,SpatialF0] = spatialFlux.compute(wvt);
        %     [StandardFp,StandardFm,StandardF0] = standardFlux.compute(wvt);
        % 
        %     self.verifyEqual(StandardFp,SpatialFp, "AbsTol",1e-7,"RelTol",1e-7);
        %     self.verifyEqual(StandardFm,SpatialFm, "AbsTol",1e-7,"RelTol",1e-7);
        %     self.verifyEqual(StandardF0,SpatialF0, "AbsTol",1e-7,"RelTol",1e-7);
        % end

        function testEnergyFluxConservation(self)
            wvt = self.wvt_;
            wvt.initWithRandomFlow(uvMax=0.1);
            
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
            self.verifyEqual(totalEnergyFlux,0, "AbsTol",1e-15);
        end

        function testTriadFluxConservation(self)
            wvt = self.wvt_;
            wvt.initWithRandomFlow(uvMax=0.1);

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

            Fp = zeros(wvt.spectralMatrixSize);
            Fm = zeros(wvt.spectralMatrixSize);
            F0 = zeros(wvt.spectralMatrixSize);
            triadFlowComponents = wvt.primaryFlowComponents;
            for i=1:length(triadFlowComponents)
                for j=1:length(triadFlowComponents)
                    [Fp_,Fm_,F0_] = wvt.nonlinearFluxForFlowComponents(triadFlowComponents(i),triadFlowComponents(j));
                    Fp = Fp + Fp_;
                    Fm = Fm + Fm_;
                    F0 = F0 + F0_;
                end
            end

            [Fp_,Fm_,F0_] = wvt.nonlinearFlux();
            self.verifyEqual(Fp,Fp_, "AbsTol",1e-15,"RelTol",1e-7);
            self.verifyEqual(Fm,Fm_, "AbsTol",1e-15,"RelTol",1e-7);
            self.verifyEqual(F0,F0_, "AbsTol",1e-15,"RelTol",1e-7);

            [Ep,Em,E0] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0,deltaT=0);
            totalEnergyFlux = sum(Ep(:))+sum(Em(:))+sum(E0(:));
            self.verifyEqual(totalEnergyFlux,0, "AbsTol",1e-15);
        end

        function testSpatialFluxConservation(self)
            wvt = self.wvt_;
            wvt.initWithRandomFlow(uvMax=0.1);

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

            wvt.addOperation(EtaTrueOperation());
            wvt.addOperation(APEOperation(wvt));
            wvt.addOperation(SpatialForcingOperation(wvt));
            int_vol = @(integrand) sum(mean(mean(shiftdim(wvt.z_int,-2).*integrand,1),2),3);

            if isa(wvt,"WVTransformHydrostatic")
                [Fu,Fv,Feta] = wvt.spatialFluxForForcingWithName("nonlinear advection");
                F_density = wvt.u .* Fu + wvt.v .* Fv+ wvt.eta_true .* shiftdim(wvt.N2,-2) .* Feta + wvt.w .* shiftdim(wvt.N2,-2) .* (wvt.eta_true-wvt.eta) ;
            elseif isa(wvt,"WVTransformBoussinesq")
                [Fu,Fv,Fw,Feta] = wvt.spatialFluxForForcingWithName("nonlinear advection");
                F_density = wvt.u .* Fu + wvt.v .* Fv +  wvt.w .* Fw + wvt.eta_true .* shiftdim(wvt.N2,-2) .* Feta + wvt.w .* shiftdim(wvt.N2,-2) .* (wvt.eta_true-wvt.eta) ;
            else
                error("Transform not yet supported.");
            end

            totalEnergyFlux = int_vol(F_density);
            self.verifyEqual(totalEnergyFlux,0, "AbsTol",1e-15);
        end

        function testNonlinearWaveTriad(self)
            % Note, this unit test is designed with the assumption that
            % Lx=Ly=2*Lz.
            wvt = self.wvt_;
            wvt.t = 0;
            L = wvt.Lx/2/pi;
            N0 = 5.2e-3;

            AbsTol = 1e-15;
            RelTol = 1e-7;

            % insert the triad
            wvt.removeAll();
            wvt.initWithWaveModes(kMode=1,lMode=0,j=1,phi=0,u=0.05,sign=1); % wave A
            wvt.addWaveModes(kMode=1,lMode=1,j=1,phi=0,u=0.05,sign=1); % wave B
            wvt.addWaveModes(kMode=0,lMode=1,j=2,phi=0,u=0.05,sign=1); % wave C

            indexA = wvt.indexFromModeNumber(1,0,1); % (2,3)
            indexB = wvt.indexFromModeNumber(1,1,1); % (2,5)
            indexC = wvt.indexFromModeNumber(0,1,2); % (3,2)

            if wvt.isHydrostatic
                Aw = sqrt(2*wvt.g/wvt.Lz)/N0;
                % equivalent depth
                ha = wvt.h_pm(2);
                hb = wvt.h_pm(2);
                hc = wvt.h_pm(3);
            else
                Aw = sqrt(2*wvt.g/((N0^2 -wvt.f^2)*wvt.Lz));

                % equivalent depth
                ha = wvt.h_pm(indexA);
                hb = wvt.h_pm(indexB);
                hc = wvt.h_pm(indexC);
            end

            N2 = N0^2/wvt.f^2;

            % non-dimensional periods of each wave
            Ta = wvt.f/wvt.Omega(indexA);
            Tb = wvt.f/wvt.Omega(indexB);
            Tc = wvt.f/wvt.Omega(indexC);

            % The true amplitude has a factor of 2 (half-complex)
            A = 2*wvt.Ap(indexA)*ha*Aw;
            B = 2*wvt.Ap(indexB)*hb*Aw/sqrt(2);
            C = 2*wvt.Ap(indexC)*hc*Aw;

            % WARNING: this unit test would be improved by testing the
            % actual forward transforms. That said, they're never actually
            % implemented in this form, so maybe this is the best we can
            % do?!
            % DCT = WVTransformConstantStratification.CosineTransformForwardMatrix(wvt.Nz);
            DCT = WVGeometryDoublyPeriodicStratifiedConstant.CosineTransformForwardMatrix(wvt.Nz);
            DCT = DCT(1:wvt.Nj,:); % dump the Nyquist mode
            % DST = WVTransformConstantStratification.SineTransformForwardMatrix(wvt.Nz);
            DST = WVGeometryDoublyPeriodicStratifiedConstant.SineTransformForwardMatrix(wvt.Nz);
            DST = cat(1,zeros(1,wvt.Nz),DST);
            DST = DST(1:wvt.Nj,:);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %% Check that our analytical expressions for (u,v,w,eta) match
            % These are pulled from the analytical expressions in my notes
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            u_bar = DCT * wvt.transformFromSpatialDomainWithFourier(wvt.u);
            v_bar = DCT * wvt.transformFromSpatialDomainWithFourier(wvt.v);
            w_bar = DST * wvt.transformFromSpatialDomainWithFourier(wvt.w);
            n_bar = DST * wvt.transformFromSpatialDomainWithFourier(wvt.eta);

            self.verifyEqual( u_bar(indexA), (A/L/2)*(-1), AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual( u_bar(indexB), (B/L/2)*(-1 + sqrt(-1)*Tb), AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual( u_bar(indexC), (C/L/2)*(-sqrt(-1)*2*Tc), AbsTol=AbsTol, RelTol=RelTol)

            self.verifyEqual( v_bar(indexA), (A/L/2)*(-sqrt(-1)*Ta), AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual( v_bar(indexB), (B/L/2)*(-1 - sqrt(-1)*Tb), AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual( v_bar(indexC), (C/L/2)*(2), AbsTol=AbsTol, RelTol=RelTol)

            self.verifyEqual( w_bar(indexA), (A/L/2)*sqrt(-1), AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual( w_bar(indexB), (B/L/2)*2*sqrt(-1), AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual( w_bar(indexC), (C/L/2)*(-sqrt(-1)), AbsTol=AbsTol, RelTol=RelTol)

            self.verifyEqual( n_bar(indexA), (A/L/2/wvt.f)*Ta, AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual( n_bar(indexB), (B/L/2/wvt.f)*2*Tb, AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual( n_bar(indexC), (C/L/2/wvt.f)*(-Tc), AbsTol=AbsTol, RelTol=RelTol)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %% Check that our analytical expressions for (uNl,vNL,wNL,etaNL) match
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            uNL = wvt.u .* wvt.diffX(wvt.u)   + wvt.v .* wvt.diffY(wvt.u)   + wvt.w .*  wvt.diffZF(wvt.u);
            vNL = wvt.u .* wvt.diffX(wvt.v)   + wvt.v .* wvt.diffY(wvt.v)   + wvt.w .*  wvt.diffZF(wvt.v);
            wNL = wvt.u .* wvt.diffX(wvt.w)   + wvt.v .* wvt.diffY(wvt.w)   + wvt.w .*  wvt.diffZG(wvt.w);
            nNL = wvt.u .* wvt.diffX(wvt.eta) + wvt.v .* wvt.diffY(wvt.eta) + wvt.w .*  wvt.diffZG(wvt.eta);
            uNL_bar = DCT * wvt.transformFromSpatialDomainWithFourier(uNL);
            vNL_bar = DCT * wvt.transformFromSpatialDomainWithFourier(vNL);
            wNL_bar = DST * wvt.transformFromSpatialDomainWithFourier(wNL);
            nNL_bar = DST * wvt.transformFromSpatialDomainWithFourier(nNL);

            % Renormalization: Must divide each amplitude by L AND because of the
            % derivative also divide by L
            uNL_A = (B*C/8/L/L/L)*(8*Tc-Tb-sqrt(-1)*(4*Tb*Tc+1));
            uNL_B = (A*C/8/L/L/L)*(-6*Tc-sqrt(-1)*(1+2*Ta*Tc));
            uNL_C = (A*B/8/L/L/L)*(Ta+Tb-sqrt(-1)*(1+Ta*Tb));

            self.verifyEqual(uNL_bar(indexA),uNL_A, AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual(uNL_bar(indexB),uNL_B, AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual(uNL_bar(indexC),uNL_C, AbsTol=AbsTol, RelTol=RelTol)

            vNL_A = (B*C/8/L/L/L)*(2*Tc - Tb + sqrt(-1)*(2*Tb*Tc-7));
            vNL_B = (A*C/8/L/L/L)*(3*Ta + sqrt(-1)*(-2*Ta*Tc-4));
            vNL_C = (A*B/8/L/L/L)*(-2*Ta-2*Tb + sqrt(-1)*(2*Ta*Tb+2));
            self.verifyEqual(vNL_bar(indexA),vNL_A, AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual(vNL_bar(indexB),vNL_B, AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual(vNL_bar(indexC),vNL_C, AbsTol=AbsTol, RelTol=RelTol)

            wNL_A = (B*C/8/L/L/L)*(5 + sqrt(-1)*(-Tb + 4*Tc));
            wNL_B = (A*C/8/L/L/L)*(-1 + sqrt(-1)*(-Ta - 2*Tc));
            wNL_C = (A*B/8/L/L/L)*(7 + sqrt(-1)*(-2*Ta-Tb));
            self.verifyEqual(wNL_bar(indexA),wNL_A, AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual(wNL_bar(indexB),wNL_B, AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual(wNL_bar(indexC),wNL_C, AbsTol=AbsTol, RelTol=RelTol)

            nNL_A = (B*C/8/L/L/L/wvt.f)*(5*Tb*Tc + sqrt(-1)*(-2*Tb + 3*Tc));
            nNL_B = (A*C/8/L/L/L/wvt.f)*(-3*Ta*Tc + sqrt(-1)*(-Ta + 2*Tc));
            nNL_C = (A*B/8/L/L/L/wvt.f)*(-Ta*Tb + sqrt(-1)*(3*Ta-4*Tb));
            self.verifyEqual(nNL_bar(indexA),nNL_A, AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual(nNL_bar(indexB),nNL_B, AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual(nNL_bar(indexC),nNL_C, AbsTol=AbsTol, RelTol=RelTol)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %% Check the equivalent dNL and zNL
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            k = shiftdim(wvt.k,-1);
            l = shiftdim(wvt.l,-1);
            dNL_bar = sqrt(-1)*(k.*uNL_bar + l.*vNL_bar);
            zNL_bar = sqrt(-1)*(k.*vNL_bar - l.*uNL_bar);

            dNL_A = (B*C/8/L/L/L/L)*(1+4*Tb*Tc+sqrt(-1)*(8*Tc-Tb));
            dNL_B = (A*C/8/L/L/L/L)*(5+4*Ta*Tc+sqrt(-1)*(3*Ta-6*Tc));
            dNL_C = (A*B/8/L/L/L/L)*(-2*Ta*Tb-2+sqrt(-1)*(-2*Ta-2*Tb));

            self.verifyEqual(dNL_bar(indexA),dNL_A, AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual(dNL_bar(indexB),dNL_B, AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual(dNL_bar(indexC),dNL_C, AbsTol=AbsTol, RelTol=RelTol)

            zNL_A = (B*C/8/L/L/L/L)*(7-2*Tb*Tc+sqrt(-1)*(2*Tc-Tb));
            zNL_B = (A*C/8/L/L/L/L)*(3+sqrt(-1)*(3*Ta+6*Tc));
            zNL_C = (A*B/8/L/L/L/L)*(-1-Ta*Tb-sqrt(-1)*(Ta+Tb));

            self.verifyEqual(zNL_bar(indexA),zNL_A, AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual(zNL_bar(indexB),zNL_B, AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual(zNL_bar(indexC),zNL_C, AbsTol=AbsTol, RelTol=RelTol)

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            %% Check the prefactors that multiply the transform
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            % signNorm = -2*(mod(wvt.j,2) == 1)+1;
            % mj = (wvt.j*pi/wvt.Lz);
            % kappa = sqrt(k.^2 + l.^2);
            % prefactor = signNorm * sqrt((wvt.g*wvt.Lz)/(2*(wvt.N0*wvt.N0 - wvt.f*wvt.f)));
            % AwD = prefactor .* (-sqrt(-1)*mj./(2*kappa));
            % AwZ = prefactor .* (-mj * wvt.f)./(2*kappa.*wvt.Omega);
            % AwW = prefactor .* (sqrt(-1)*kappa./2);
            % AwN = prefactor .* (-(wvt.N0*wvt.N0)*kappa./(2*wvt.Omega));

            % Fp = AwD .* dNL_bar + AwZ .* zNL_bar + AwW .* wNL_bar + AwN .* nNL_bar;
            % Ep = 2*wvt.Apm_TE_factor.*real( Fp .* conj(wvt.Ap) );
            [Fp,Fm,F0] = wvt.nonlinearFlux();
            [Ep,~,~] = wvt.energyFluxFromNonlinearFlux(Fp,Fm,F0,deltaT=0);

            if wvt.isHydrostatic
                E_A = -pi*(A*B*C/16/L/L/L)*(7*Ta +Tb - 8*Tc + (5*N2 - 2)*Ta*Tb*Tc);
                E_B = -pi*(A*B*C/16/L/L/L)*(-3*Ta + 3*Tb + 6*Tc - 6*N2*Ta*Tb*Tc);
                E_C = -pi*(A*B*C/16/L/L/L)*(-4*Ta - 4*Tb + 2*Tc  +(N2+2)*Ta*Tb*Tc);
            else
                E_A = -pi*(A*B*C/16/L/L/L)*(7*Ta - 4*Tc + (5*N2 - 2)*Ta*Tb*Tc);
                E_B = -pi*(A*B*C/16/L/L/L)*(-5*Ta +3*Tb + 2*Tc - 6*N2*Ta*Tb*Tc);
                E_C = -pi*(A*B*C/16/L/L/L)*(-2*Ta - 3*Tb + 2*Tc  +(N2+2)*Ta*Tb*Tc);
            end
            self.verifyEqual(Ep(indexA),E_A, AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual(Ep(indexB),E_B, AbsTol=AbsTol, RelTol=RelTol)
            self.verifyEqual(Ep(indexC),E_C, AbsTol=AbsTol, RelTol=RelTol)

        end

    end

end