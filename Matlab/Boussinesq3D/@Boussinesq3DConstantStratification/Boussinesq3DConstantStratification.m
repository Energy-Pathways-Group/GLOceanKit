classdef Boussinesq3DConstantStratification < handle
    %3D Boussinesq model with constant stratification solved in wave-vortex
    %space
    
    properties
        x, y, z
        k, l, j
        Nx, Ny, Nz
        Lx, Ly, Lz
        f0, N0, rho0
        iOmega
        
        realScratch, complexScratch; % of size Nx x Ny x (2*Nz-1)
        F,G
        ApU, ApV, ApN
        AmU, AmV, AmN
        A0U, A0V, A0N
        
        UAp, UAm, UA0
        VAp, VAm, VA0
        WAp, WAm
        NAp, NAm, NA0
        
        dt
        integrator
        t0 = 0 % current time that the coefficients {Ap,Am,A0} are wound to
        t = 0 % current time of the integrator
        Y % cell array with {Ap,Am,A0}
        
        nParticles = 0
        shouldAntialias = 0;
    end
    
    properties (Dependent)
        % These convert the coefficients to their depth integrated energies
        Apm_HKE_factor
        Apm_VKE_factor
        Apm_PE_factor
        Apm_TE_factor

        A0_HKE_factor
        A0_PE_factor
        A0_TE_factor
    end
    
    properties (Constant)
        g = 9.81;
    end
    
    methods
        function self = Boussinesq3DConstantStratification(dims, n, latitude, N0, rho0)
            % rho0 is optional.
            if length(dims) ~=3 || length(n) ~= 3
                error('The dims and n variables must be of length 3. You need to specify x,y,z');
            end
            
            if mod(log2(n(3)),1) == 0
                error('You are implicitly asking for periodic boundary conditions in the vertical. This is not supported.');
            elseif mod(log2(n(3)-1),1) == 0
                self.Nz = n(3);
            else
                error('The vertical dimension must have 2^n or (2^n)+1 points. This is an artificial restriction.');
            end
                        
            self.Lx = dims(1);
            self.Ly = dims(2);
            self.Lz = dims(3);
            
            self.Nx = n(1);
            self.Ny = n(2);
            self.Nz = n(3);

            dx = self.Lx/self.Nx;
            dy = self.Ly/self.Ny;
            dz = self.Lz/(self.Nz-1);
            
            self.x = dx*(0:self.Nx-1)'; % periodic basis
            self.y = dy*(0:self.Ny-1)'; % periodic basis
            self.z = dz*(0:(self.Nz-1))' - self.Lz; % Cosine basis for DCT-I and DST-I
                        
            dk = 1/self.Lx;          % fourier frequency
            self.k = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]*dk)';
            dl = 1/self.Ly;          % fourier frequency
            self.l = 2*pi*([0:ceil(self.Ny/2)-1 -floor(self.Ny/2):-1]*dl)';
            self.j = (0:(self.Nz-2))';
            
            self.N0 = N0;
            self.f0 = 2 * 7.2921E-5 * sin( latitude*pi/180 );
            if ~exist('rho0','var')
                self.rho0 = 1025;
            else
                self.rho0 = rho0;
            end
            
            self.BuildTransformationMatrices();
            
            % Preallocate this array for a faster dct
            self.realScratch = zeros(self.Nx,self.Ny,(2*self.Nz-1));
            self.complexScratch = complex(zeros(self.Nx,self.Ny,2*(self.Nz-1)));
            
            % Now set the initial conditions to zero
            Ap0 = zeros(size(self.ApU));
            Am0 = zeros(size(self.ApU));
            A00 = zeros(size(self.ApU));
            self.Y = {Ap0;Am0;A00;};
            
%             self.integrator = ArrayIntegrator(@(t,y0) self.NonlinearFluxAtTimeArray(t,y0),self.Y,2*pi/self.f0/100);
        end
        
        function IncrementForward(self)
            self.integrator.currentY = self.Y;
            self.Y = self.integrator.IncrementForward();
            self.t = self.integrator.currentTime;
        end
        
        function h = h(self)
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            K2 = K.*K + L.*L;
            M = J*pi/self.Lz;
            h = (1/self.g)*(self.N0*self.N0 - self.f0*self.f0)./(M.*M+K2);
            h(:,:,1) = 1; % prevent divide by zero
        end
        
        function Kh = Kh(self)
            [K,L,~] = ndgrid(self.k,self.l,self.j);
            Kh = sqrt(K.*K + L.*L);
        end
        
        function Omega = Omega(self)
            [K,L,~] = ndgrid(self.k,self.l,self.j);
            Omega = sqrt(self.g*self.h.*(K.*K + L.*L) + self.f0*self.f0);
        end
        
        function self = BuildTransformationMatrices(self)
            % Build wavenumbers
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            alpha = atan2(L,K);
            K2 = K.*K + L.*L;
            Kh = sqrt(K2);      % Total horizontal wavenumber
            M = J*pi/self.Lz;        % Vertical wavenumber
            
            N = self.N0;
            f = self.f0;
            
            g_ = 9.81;
            h = (1/g_)*(N*N-f*f)./(M.*M+K2);
            h(:,:,1) = 1; % prevent divide by zero
            
            omega = sqrt(g_*h.*K2 + f*f);
            if abs(self.f0) < 1e-14 % This handles the f=0 case.
                omega(omega == 0) = 1;
            end
            fOmega = f./omega;
            self.iOmega = InternalWaveModel.MakeHermitian(sqrt(-1)*omega);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Normalization for the vertical modes
            % This comes from equations B12 in the manuscript.
            signNorm = -2*(mod(J,2) == 1)+1; % equivalent to (-1)^j
            self.F = signNorm .* (h.*M)*sqrt(2*g_/(self.Lz*(N*N-f*f)));
            self.G = signNorm .* sqrt(2*g_/(self.Lz*(N*N-f*f)));
            self.F(:,:,1) = 2; % j=0 mode is a factor of 2 too big in DCT-I
            self.G(:,:,1) = 1; % j=0 mode doesn't exist for G
            
            MakeHermitian = @(f) InternalWaveModel.MakeHermitian(f);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Transform matrices (U,V,N) -> (Ap,Am,A0)
            % This comes from equations B13 and B14 in the manuscript
            % or equation 5.5 without the factor of h.
            self.ApU = (1/2)*(cos(alpha)+sqrt(-1)*fOmega.*sin(alpha));
            self.ApV = (1/2)*(sin(alpha)-sqrt(-1)*fOmega.*cos(alpha));
            self.ApN = -g_*Kh./(2*omega);
            
            self.AmU = (1/2)*(cos(alpha)-sqrt(-1)*fOmega.*sin(alpha));
            self.AmV = (1/2)*(sin(alpha)+sqrt(-1)*fOmega.*cos(alpha));
            self.AmN = g_*Kh./(2*omega);
            
            % There are no k^2+l^2>0, j=0 wave solutions. Only the inertial
            % solution exists at k=l=j=0.
            self.ApU(:,:,1) = 0;
            self.ApV(:,:,1) = 0;
            self.ApN(:,:,1) = 0;
            
            self.AmU(:,:,1) = 0;
            self.AmV(:,:,1) = 0;
            self.AmN(:,:,1) = 0;
            
            % Now set the inertial stuff (this is just a limit of above)
            self.ApU(1,1,:) = 1/2;
            self.ApV(1,1,:) = -sqrt(-1)/2;
            self.AmU(1,1,:) = 1/2;
            self.AmV(1,1,:) = sqrt(-1)/2;
            
            % Equation B14
            self.A0U = sqrt(-1)*h.*(fOmega./omega) .* L;
            self.A0V = -sqrt(-1)*h.*(fOmega./omega) .* K;
            self.A0N = fOmega.^2;
            
            % k > 0, l > 0, j=0; Equation B11 in the manuscript
            self.A0U(:,:,1) =  sqrt(-1)*(f/g_)*L(:,:,1)./K2(:,:,1); % Note the divide by zero at k=l=0
            self.A0V(:,:,1) = -sqrt(-1)*(f/g_)*K(:,:,1)./K2(:,:,1);
            self.A0N(:,:,1) = 0;
            
            % The k=l=0, j>=0 geostrophic solutions are a simple density anomaly
            self.A0U(1,1,:) = 0;
            self.A0V(1,1,:) = 0;
            self.A0N(1,1,:) = 1;
            self.A0N(1,1,1) = 0;
            
            % Finally, we need to take care of the extra factor of 2 that
            % comes out of the discrete cosine transform
            
            % http://helper.ipam.ucla.edu/publications/mtws1/mtws1_12187.pdf
            shouldAntiAlias = 1;
            AntiAliasFilter = ones(size(self.ApU));
            if shouldAntiAlias == 1
                AntiAliasFilter(Kh > 2*max(abs(self.k))/3 | J > 2*max(abs(self.j))/3) = 0;
            end
            
            % Now make the Hermitian conjugate match.
            iFTransformScaling = 2./(self.Nx*self.Ny*self.F);
            iGTransformScaling = 2./(self.Nx*self.Ny*self.G);
            self.ApU = AntiAliasFilter .* iFTransformScaling .* MakeHermitian(self.ApU);
            self.ApV = AntiAliasFilter .* iFTransformScaling .* MakeHermitian(self.ApV);
            self.ApN = AntiAliasFilter .* iGTransformScaling .* MakeHermitian(self.ApN);
            
            self.AmU = AntiAliasFilter .* iFTransformScaling .* MakeHermitian(self.AmU);
            self.AmV = AntiAliasFilter .* iFTransformScaling .* MakeHermitian(self.AmV);
            self.AmN = AntiAliasFilter .* iGTransformScaling .* MakeHermitian(self.AmN);
            
            self.A0U = AntiAliasFilter .* iFTransformScaling .* MakeHermitian(self.A0U);
            self.A0V = AntiAliasFilter .* iFTransformScaling .* MakeHermitian(self.A0V);
            self.A0N = AntiAliasFilter .* iGTransformScaling .* MakeHermitian(self.A0N);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Transform matrices (Ap,Am,A0) -> (U,V,W,N)
            % These can be pulled from equation C4 in the manuscript
            self.UAp = (cos(alpha)-sqrt(-1)*fOmega.*sin(alpha));
            self.UAm = (cos(alpha)+sqrt(-1)*fOmega.*sin(alpha));
            self.UA0 = -sqrt(-1)*(g_/f)*L;
            
            self.VAp = (sin(alpha)+sqrt(-1)*fOmega.*cos(alpha));
            self.VAm = (sin(alpha)-sqrt(-1)*fOmega.*cos(alpha));
            self.VA0 = sqrt(-1)*(g_/f)*K;
                
            self.WAp = -sqrt(-1)*Kh.*h;
            self.WAm = -sqrt(-1)*Kh.*h;
            
            self.NAp = -Kh.*h./omega;
            self.NAm = Kh.*h./omega;
            self.NA0 = ones(size(Kh));
            
            % No buoyancy anomaly for j=0 geostrophic solutions
            self.NA0(:,:,1) = 0;
            
            % There are no k^2+l^2>0, j=0 wave solutions. 
            self.UAp(:,:,1) = 0;
            self.VAp(:,:,1) = 0;
            self.NAp(:,:,1) = 0;
            
            self.UAm(:,:,1) = 0;
            self.VAm(:,:,1) = 0;
            self.NAm(:,:,1) = 0;
            
            % Only the inertial solution exists at k=l=j=0 as a negative
            % wave.
            self.UAp(1,1,:) = 1;
            self.VAp(1,1,:) = sqrt(-1);
            self.UAm(1,1,:) = 1;
            self.VAm(1,1,:) = -sqrt(-1);
            
            if abs(self.f0) < 1e-14 % This handles the f=0 case.
                self.UA0 = zeros(size(Kh));
                self.VA0 = zeros(size(Kh));
                self.NA0 = zeros(size(Kh));
            end
            
            % Now make the Hermitian conjugate match AND pre-multiply the
            % coefficients for the transformations.
            FTransformScaling = 0.5*self.Nx*self.Ny*self.F;
            self.UAp = FTransformScaling .* MakeHermitian(self.UAp);
            self.UAm = FTransformScaling .* MakeHermitian(self.UAm);
            self.UA0 = FTransformScaling .* MakeHermitian(self.UA0);
            
            self.VAp = FTransformScaling .* MakeHermitian(self.VAp);
            self.VAm = FTransformScaling .* MakeHermitian(self.VAm);
            self.VA0 = FTransformScaling .* MakeHermitian(self.VA0);
            
            GTransformScaling = 0.5*self.Nx*self.Ny*self.G;
            self.WAp = GTransformScaling .* MakeHermitian(self.WAp);
            self.WAm = GTransformScaling .* MakeHermitian(self.WAm);
            
            self.NAp = GTransformScaling .* MakeHermitian(self.NAp);
            self.NAm = GTransformScaling .* MakeHermitian(self.NAm);
            self.NA0 = GTransformScaling .* MakeHermitian(self.NA0);
        end
          
        function [Ap,Am,A0] = TransformUVEtaToWaveVortex(self,U,V,N)
            % This is the 'S^{-1}' operator (C5) in the manuscript
            Ubar = self.TransformFromSpatialDomainWithF(U);
            Vbar = self.TransformFromSpatialDomainWithF(V);
            Nbar = self.TransformFromSpatialDomainWithG(N);
            
            Ap = self.ApU.*Ubar + self.ApV.*Vbar + self.ApN.*Nbar;
            Am = self.AmU.*Ubar + self.AmV.*Vbar + self.AmN.*Nbar;
            A0 = self.A0U.*Ubar + self.A0V.*Vbar + self.A0N.*Nbar;
        end
        
        function [U,V,W,N] = TransformWaveVortexToUVWEta(self,Ap,Am,A0)
            % This is the 'S' operator (C4) in the manuscript
            Ubar = self.UAp.*Ap + self.UAm.*Am + self.UA0.*A0;
            Vbar = self.VAp.*Ap + self.VAm.*Am + self.VA0.*A0;
            Wbar = self.WAp.*Ap + self.WAm.*Am;
            Nbar = self.NAp.*Ap + self.NAm.*Am + self.NA0.*A0;
            
            U = self.TransformToSpatialDomainWithF(Ubar);
            V = self.TransformToSpatialDomainWithF(Vbar);
            W = self.TransformToSpatialDomainWithG(Wbar);
            N = self.TransformToSpatialDomainWithG(Nbar);
        end
        
        function [uNL,vNL,nNL] = NonlinearFluxFromSpatial(self,u,v,w,eta)
            uNL = u.*DiffFourier(self.x,u,1,1) + v.*DiffFourier(self.y,u,1,2) + w.*DiffCosine(self.z,u,1,3);
            vNL = u.*DiffFourier(self.x,v,1,1) + v.*DiffFourier(self.y,v,1,2) + w.*DiffCosine(self.z,v,1,3);
            nNL = u.*DiffFourier(self.x,eta,1,1) + v.*DiffFourier(self.y,eta,1,2) + w.*DiffSine(self.z,eta,1,3);
        end
        
        function F = NonlinearFluxAtTimeArray(self,t,Y0)
            [Fp,Fm,F0] = self.NonlinearFluxAtTime(t,Y0{1},Y0{2},Y0{3});
            F = {Fp,Fm,F0};
        end
        
        function [Fp,Fm,F0] = NonlinearFluxAtTime(self,t,Ap,Am,A0)
            % Apply operator T_\omega---defined in (C2) in the manuscript
            phase = exp(self.iOmega*(t-self.t0));
            Ap = Ap .* phase;
            Am = Am .* conj(phase);
            
            % Apply operator S---defined in (C4) in the manuscript
            Ubar = self.UAp.*Ap + self.UAm.*Am + self.UA0.*A0;
            Vbar = self.VAp.*Ap + self.VAm.*Am + self.VA0.*A0;
            Wbar = self.WAp.*Ap + self.WAm.*Am;
            Nbar = self.NAp.*Ap + self.NAm.*Am + self.NA0.*A0;
            
            % Finishing applying S, but also compute derivatives at the
            % same time
            [U,Ux,Uy,Uz] = self.TransformToSpatialDomainWithFAllDerivatives(Ubar);
            [V,Vx,Vy,Vz] = self.TransformToSpatialDomainWithFAllDerivatives(Vbar);
            W = self.TransformToSpatialDomainWithG(Wbar);
            [ETA,ETAx,ETAy,ETAz] = self.TransformToSpatialDomainWithGAllDerivatives(Nbar);
            
            % Compute the nonlinear terms in the spatial domain
            % (pseudospectral!)
            uNL = -U.*Ux - V.*Uy - W.*Uz;
            vNL = -U.*Vx - V.*Vy - W.*Vz;
            nNL = -U.*ETAx - V.*ETAy - W.*ETAz;
            
            % Now apply the operator S^{-1} and then T_\omega^{-1}
            uNLbar = self.TransformFromSpatialDomainWithF(uNL);
            vNLbar = self.TransformFromSpatialDomainWithF(vNL);
            nNLbar = self.TransformFromSpatialDomainWithG(nNL);
            
%             CheckHermitian(uNLbar+vNLbar);
%             CheckHermitian(nNLbar);

            Fp = (self.ApU.*uNLbar + self.ApV.*vNLbar + self.ApN.*nNLbar) .* conj(phase);
            Fm = (self.AmU.*uNLbar + self.AmV.*vNLbar + self.AmN.*nNLbar) .* phase;
            F0 = self.A0U.*uNLbar + self.A0V.*vNLbar + self.A0N.*nNLbar;
            
%             CheckHermitian(Fp+Fm);
%             CheckHermitian(F0);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Energetics
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function value = get.Apm_TE_factor(self)
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            K2 = K.*K + L.*L;
            M = J*pi/self.Lz;
            h = (1/self.g)*(self.N0*self.N0-self.f0*self.f0)./(M.*M+K2);
            
            value = h; % factor of 2 larger than in the manuscript
            value(:,:,1) = self.Lz;
        end
        
        function value = get.A0_HKE_factor(self)
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            K2 = K.*K + L.*L;
            M = J*pi/self.Lz;
            h = (1/self.g)*(self.N0*self.N0-self.f0*self.f0)./(M.*M+K2);
            h(:,:,1) = 1; % prevent divide by zero 
            
            % This comes from equation (3.10) in the manuscript, but using
            % the relation from equation A2b
            % omega = sqrt(self.g*h.*K2 + self.f0*self.f0);
            % value = (self.g/(self.f0*self.f0)) * (omega.*omega - self.f0*self.f0) .* (self.N0*self.N0 - omega.*omega) / (2 * (self.N0*self.N0 - self.f0*self.f0) );
            value = (self.g^3/(self.f0*self.f0)) * K2.*h.*h.*M.*M / (2 * (self.N0*self.N0 - self.f0*self.f0) ); % factor of 2 larger than in the manuscript
            value(:,:,1) = (self.g^2/(self.f0*self.f0)) * K2(:,:,1) * self.Lz/2;
        end
        function value = get.A0_PE_factor(self)
            value = self.g*self.N0*self.N0/(self.N0*self.N0-self.f0*self.f0)/2; % factor of 2 larger than in the manuscript
        end
        function value = get.A0_TE_factor(self)
            value = self.A0_HKE_factor + self.A0_PE_factor;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Energetics (total)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function energy = totalEnergy(self)
            [u,v,w,eta] = self.VariableFieldsAtTime(0,'u','v','w','eta');
            energy = trapz(self.z,mean(mean( u.^2 + v.^2 + w.^2 + self.N0*self.N0*eta.*eta, 1 ),2 ) )/2;
        end
        
        function energy = totalSpectralEnergy(self)
%             energy = self.inertialEnergy + self.waveEnergy + self.geostrophicEnergy;
App = self.Ap; Amm = self.Am; A00 = self.A0;
energy = sum(sum(sum( self.Apm_TE_factor.*( App.*conj(App) + Amm.*conj(Amm) ) + self.A0_TE_factor.*( A00.*conj(A00) ) )));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Major constituents
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function energy = inertialEnergy(self)
            energy = self.barotropicInertialEnergy + self.baroclinicInertialEnergy;
        end
        
        function energy = waveEnergy(self)
            energy = self.internalWaveEnergyPlus + self.internalWaveEnergyMinus;
        end
        
        function energy = geostrophicEnergy(self)
            energy = self.barotropicGeostrophicEnergy + self.baroclinicGeostrophicEnergy;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Geostrophic constituents
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function energy = barotropicGeostrophicEnergy(self)
            C = self.A0_TE_factor;
            B = self.A0;
            energy = sum(sum(sum( C(:,:,1) .* (B(:,:,1).*conj(B(:,:,1))) )));
        end
        
        function energy = baroclinicGeostrophicEnergy(self)
            C = self.A0_TE_factor;
            B = self.A0;
            energy = sum(sum(sum( C(:,:,2:end) .* (B(:,:,2:end).*conj(B(:,:,2:end))) )));
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Inertia-gravity wave constituents
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function energy = barotropicInertialEnergy(self)
            App = self.Ap;
            Amm = self.Am;
            C = self.Apm_TE_factor;
            energy = C(1,1,1)*( App(1,1,1).*conj(App(1,1,1)) + Amm(1,1,1).*conj(Amm(1,1,1)) );
        end
        
        function energy = baroclinicInertialEnergy(self)
            App = self.Ap;
            Amm = self.Am;
            C = self.Apm_TE_factor;
            energy = sum(sum(sum( C(1,1,2:end).* (abs(App(1,1,2:end)).^2 + abs(Amm(1,1,2:end)).^2) )));
        end
        
        function energy = internalWaveEnergyPlus(self)
            A = self.Ap;
            A(1,1,:) = 0;
            C = self.Apm_TE_factor;
            energy = sum( C(:).* (A(:).*conj(A(:)))  );
        end
        
        function energy = internalWaveEnergyMinus(self)
            A = self.Am;
            A(1,1,:) = 0;
            C = self.Apm_TE_factor;
            energy = sum( C(:).* (A(:).*conj(A(:)))  );
        end
        
        function summarizeEnergyContent(self)
            total = self.totalSpectralEnergy;
            ioPct = 100*self.inertialEnergy/total;
            wavePct = 100*self.waveEnergy/total;
            gPct = 100*self.geostrophicEnergy/total;
            wavePlusPct = 100*self.internalWaveEnergyPlus/self.waveEnergy;
            waveMinusPct = 100*self.internalWaveEnergyMinus/self.waveEnergy;
            
            fprintf('%.1f m^3/s^2 total depth integrated energy, split (%.1f,%.1f,%.1f) between (inertial,wave,geostrophic) with wave energy split %.1f/%.1f +/-\n',total,ioPct,wavePct,gPct,wavePlusPct,waveMinusPct);
        end
        
        
        function Ap = Ap(self)
            Ap = self.Y{1};
        end
        function Am = Am(self)
            Am = self.Y{2};
        end
        function A0 = A0(self)
            A0 = self.Y{3};
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Create a single wave (public)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        period = InitializeWithPlaneWave(self, k0, l0, j0, UAmp, sign)
        
        RemoveAllGriddedWaves(self)
        
        [omega,k,l] = SetGriddedWavesWithWavemodes(self, kMode, lMode, jMode, phi, Amp, signs)
        
        [omega,k,l] = AddGriddedWavesWithWavemodes(self, kMode, lMode, jMode, phi, Amp, signs)
           
        
        function ratio = UmaxGNormRatioForWave(self,k0, l0, j0)
            if j0 == 0
                ratio = 1;
            else
                ratio = abs(1/self.F(k0+1,l0+1,j0+1));
            end
        end     
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Validation and internal unit testing
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [C11,C21,C31,C12,C22,C32,C13,C23,C33] = Validate(self)
            % This is S*S^{-1} and therefore returns the values in
            % wave-vortex space. So, C11 represents Ap and should be 1s
            % where we expected Ap solutions to exist.
            C11 = self.ApU.*self.UAp + self.ApV.*self.VAp + self.ApN.*self.NAp;
            C21 = self.AmU.*self.UAp + self.AmV.*self.VAp + self.AmN.*self.NAp;
            C31 = self.A0U.*self.UAp + self.A0V.*self.VAp + self.A0N.*self.NAp;
            
            C12 = self.ApU.*self.UAm + self.ApV.*self.VAm + self.ApN.*self.NAm;
            C22 = self.AmU.*self.UAm + self.AmV.*self.VAm + self.AmN.*self.NAm;
            C32 = self.A0U.*self.UAm + self.A0V.*self.VAm + self.A0N.*self.NAm;
            
            C13 = self.ApU.*self.UA0 + self.ApV.*self.VA0 + self.ApN.*self.NA0;
            C23 = self.AmU.*self.UA0 + self.AmV.*self.VA0 + self.AmN.*self.NA0;
            C33 = self.A0U.*self.UA0 + self.A0V.*self.VA0 + self.A0N.*self.NA0;
            
            % Maybe check that max(abs(C12(:))) is very small.
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations to and from the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function u_bar = TransformFromSpatialDomainWithF(self, u)
            % All coefficients are subsumbed into the transform
            % coefficients ApU,ApV,etc.
%             self.dctScratch = ifft(cat(3,u,u(:,:,(self.Nz-1):-1:2)),2*(self.Nz-1),3);
%             u_bar = real(self.dctScratch(:,:,1:(self.Nz-1))); % include barotropic mode, but leave off the Nyquist.
%             u_bar = fft(fft(u_bar,self.Nx,1),self.Ny,2);
%             
            u = ifft(cat(3,u,flip(u,3)),2*(self.Nz-1),3,'symmetric');
            u_bar = fft(fft(u(:,:,1:(self.Nz-1)),self.Nx,1),self.Ny,2);
        end
        
        function w_bar = TransformFromSpatialDomainWithG(self, w)
            % df = 1/(2*(Nz-1)*dz)
            % nyquist = (Nz-2)*df
            w = ifft(cat(3,w,-w(:,:,(self.Nz-1):-1:2)),2*(self.Nz-1),3);
            w_bar = imag(w(:,:,1:(self.Nz-1)));
            w_bar = fft(fft(w_bar,self.Nx,1),self.Ny,2);
        end
        
        function u = TransformToSpatialDomainWithF(self, u_bar)
            % All coefficients are subsumbed into the transform
            % coefficients UAp,UAm,etc.
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');    
            
            % Re-order to convert to a DCT-I via FFT
            u = cat(3,u,zeros(self.Nx,self.Ny),flip(u,3));
            u = ifft( u,2*(self.Nz-1),3,'symmetric');
            u = u(:,:,1:self.Nz)*(2*self.Nz-2); % We do not incorporate this coefficient into UAp, etc, so that the transforms remain inverses
        end  
                
        function w = TransformToSpatialDomainWithG(self, w_bar )
            % All coefficients are subsumbed into the transform
            % coefficients NAp,NAm,etc.
            w = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');
            
            % Re-order to convert to a DST-I via FFT
            w = sqrt(-1)*cat(3,-w,zeros(self.Nx,self.Ny),flip(w,3));
            w = ifft( w,2*(self.Nz-1),3,'symmetric');
            w = w(:,:,1:self.Nz)*(2*self.Nz-2); % We do not incorporate this coefficient into UAp, etc, so that the transforms remain inverses
        end
        
        function [u,ux,uy,uz] = TransformToSpatialDomainWithFAllDerivatives(self, u_bar)
            % All coefficients are subsumbed into the transform
            % coefficients UAp,UAm,etc.
            u = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');    
            
            % Re-order to convert to a DCT-I via FFT
            self.realScratch = cat(3,u,zeros(self.Nx,self.Ny),flip(u,3));
            u = ifft( self.realScratch,2*(self.Nz-1),3,'symmetric');
            u = u(:,:,1:self.Nz)*(2*self.Nz-2);
            
            ux = ifft( sqrt(-1)*self.k.*fft(u,self.Nx,1), self.Nx, 1,'symmetric');
            uy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(u,self.Ny,2), self.Ny, 2,'symmetric');
            
            % To take the derivative, we multiply, then sine transform back
            m = reshape(-pi*self.j/self.Lz,1,1,[]);
            % That would be the normal multiplier, but we want to get it
            % ready for a DST, so.. -i*m,0,i*flip(m)
            m = cat(3,-sqrt(-1)*m,0,sqrt(-1)*flip(m,3));
            uz = ifft( m.*self.realScratch,2*(self.Nz-1),3,'symmetric');
            uz = uz(:,:,1:self.Nz)*(2*self.Nz-2);
        end  
        
        function [w,wx,wy,wz] = TransformToSpatialDomainWithGAllDerivatives(self, w_bar )
            % All coefficients are subsumbed into the transform
            % coefficients NAp,NAm,etc.
            w = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');
            
            % Re-order to convert to a DST-I via FFT
            self.complexScratch = sqrt(-1)*cat(3,-w,zeros(self.Nx,self.Ny),flip(w,3));
            w = ifft( self.complexScratch,2*(self.Nz-1),3,'symmetric');
            w = w(:,:,1:self.Nz)*(2*self.Nz-2); % We do not incorporate this coefficient into UAp, etc, so that the transforms remain inverses
            
            wx = ifft( sqrt(-1)*self.k.*fft(w,self.Nx,1), self.Nx, 1,'symmetric');
            wy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(w,self.Ny,2), self.Ny, 2,'symmetric');
            
            % To take the derivative, we multiply, then cosine transform back
            m = reshape(pi*self.j/self.Lz,1,1,[]);
            % That would be the normal multiplier, but we want to get it
            % ready for a DST, so.. i*m,0,-i*flip(m)
            m = cat(3,sqrt(-1)*m,0,-sqrt(-1)*flip(m,3));
            wz = ifft( m.*self.complexScratch,2*(self.Nz-1),3,'symmetric');
            wz = wz(:,:,1:self.Nz)*(2*self.Nz-2);
        end
        
        function [ApIO,AmIO,ApIGW,AmIGW,A0G,A0G0,A0rhobar] = GenerateRandomFlowState(self)
            % Random flow state, separated out by solution type.
            % Adding the solution types together, gives a complete state.
            % Ap = ApIO + ApIGW;
            % Am = AmIO + AmIGW;
            % A0 = A0G + A0G0 + A0rhobar;
            shouldExcludeNyquist = 1;
            ApIGW = Boussinesq3DConstantStratification.GenerateHermitianRandomMatrix( size(self.G), shouldExcludeNyquist );
            AmIGW = Boussinesq3DConstantStratification.GenerateHermitianRandomMatrix( size(self.G), shouldExcludeNyquist );
            A0G = 6e-2*Boussinesq3DConstantStratification.GenerateHermitianRandomMatrix( size(self.G), shouldExcludeNyquist );
            
            ApIO = zeros(size(self.G));
            AmIO = zeros(size(self.G));
            A0G0 = zeros(size(self.G));
            A0rhobar = zeros(size(self.G));
            
            % inertial oscillations only exist at k=l=0
            ApIO(1,1,:) = ApIGW(1,1,:);
            AmIO(1,1,:) = conj(ApIGW(1,1,:));
            
            % zero out all j=0, and k=l=0 values.
            ApIGW(:,:,1) = 0;
            ApIGW(1,1,:) = 0;
            AmIGW(:,:,1) = 0;
            AmIGW(1,1,:) = 0;
            
            % barotropic geostrophic at all k and l>0, j=0
            A0G0(:,:,1) = 0.1*A0G(:,:,1);
            A0G0(1,1,1) = 0;
            
            % mean density anomaly
            A0rhobar(1,1,2:end) = real(A0G(1,1,2:end));
            
            % zero out all j=0, and k=l=0 values.
            A0G(1,1,:) = 0;
            A0G(:,:,1) = 0;
        end
        
    end
    
    methods (Static)
        function A = CheckHermitian(A)
            M = size(A,1);
            N = size(A,2);
            K = size(A,3);
            
            for k=1:K
                for i=M:-1:1
                    for j=N:-1:1
                        ii = mod(M-i+1, M) + 1;
                        jj = mod(N-j+1, N) + 1;
                        if A(i,j,k) ~= conj(A(ii,jj,k))
                            error('(i,j,k)=(%d,%d,%d) is not conjugate with (%d,%d,%d)\n',i,j,k,ii,jj,k)
                        end
                    end
                end
            end
            
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Forces a 3D matrix to be Hermitian, ready for an FFT (internal)
        %
        % The approach taken here is that the (k=-Nx/2..Nx/2,l=0..Ny/2+1) wave
        % numbers are primary, and the (k=-Nx/2..Nx/2,l=-Ny/2..1) are inferred as
        % conjugates. Also, the negative k wavenumbers for l=0. The Nyquist wave
        % numbers are set to zero to avoid complications.
        %
        % This function is NOT a true "Make Hermitian" function because it
        % doesn't force the k=l=0 to be real.
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function A = MakeHermitian(A)
            M = size(A,1);
            N = size(A,2);
            K = size(A,3);
            
            % The order of the for-loop is chosen carefully here.
            for k=1:K
                for j=1:(N/2+1)
                    for i=1:M
                        ii = mod(M-i+1, M) + 1;
                        jj = mod(N-j+1, N) + 1;
                        if i == ii && j == jj
                            % A(i,j,k) = real(A(i,j,k)); % self-conjugate term
                            % This is not normally what you'd do, but we're being
                            % tricky by later adding the conjugate
                            if i == 1 && j == 1
                                continue;
                            else
                                A(i,j,k) = 0;
                            end
                        elseif j == N/2+1 % Kill the Nyquist, rather than fix it.
                            A(i,j,k) = 0;
                        else % we are letting l=0, k=Nx/2+1 terms set themselves again, but that's okay
                            A(ii,jj,k) = conj(A(i,j,k));
                        end
                    end
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Returns a matrix the same size as A with 1s at the 'redundant'
        % hermiation indices.
        function A = RedundantHermitianCoefficients(A)
            A = zeros(size(A));
            
            M = size(A,1);
            N = size(A,2);
            K = size(A,3);
            
            % The order of the for-loop is chosen carefully here.
            for k=1:K
                for j=1:(N/2+1)
                    for i=1:M
                        ii = mod(M-i+1, M) + 1;
                        jj = mod(N-j+1, N) + 1;
                        if i == ii && j == jj
                            if i == 1 && j == 1
                                continue;
                            else
                                A(i,j,k) = 0;
                            end
                        elseif j == N/2+1
                            A(i,j,k) = 0;
                        else
                            if j == 1 && i > M/2
                                continue;
                            end
                            A(ii,jj,k) = 1;
                        end
                    end
                end
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Returns a matrix the same size as A with 1s at the Nyquist
        % frequencies.
        function A = NyquistWavenumbers(A)
            A = zeros(size(A));
            Nx = size(A,1);
            Ny = size(A,2);
            if Nx > 1
                A(ceil(Nx/2)+1,:) = 1;
            end
            if Ny > 1
                A(:,ceil(Ny/2)+1) = 1;
            end
        end
        
        function A = GenerateHermitianRandomMatrix( size, shouldExcludeNyquist )
            
            nX = size(1); nY = size(2);
            if length(size) > 2
                nZ = size(3);
            else
                nZ = 1;
            end
            A = InternalWaveModel.MakeHermitian(randn(size) + sqrt(-1)*randn(size) )/sqrt(2);
            if shouldExcludeNyquist == 1
                mask = ~InternalWaveModel.NyquistWavenumbers(A(:,:,1));
                A = mask.*A;
            else
                A(1,1,:) = 2*A(1,1,:); % Double the zero frequency
                A(nX/2+1,1,:) = -2*real(A(nX/2+1,1,:)); % Double the Nyquist frequency
                A(1,nY/2+1,:) = -2*real(A(1,nY/2+1,:)); % Double the Nyquist frequency
                A(nX/2+1,nY/2+1,:) = -2*real(A(nX/2+1,nY/2+1,:)); % Double the Nyquist frequency
            end
            
        end
    end
        
        
end 



