classdef WaveVortexModel < handle
    %3D Boussinesq model with constant stratification solved in wave-vortex
    %space
    
    properties
        x, y, z
        k, l, j

        X,Y,Z

        Nx, Ny, Nz, nModes
        Nk, Nl, Nj % actual sizes in the spectral domain
        Lx, Ly, Lz
        f0, Nmax, rho0, latitude
        iOmega
        rhobar, N2, dLnN2 % on the z-grid, size(N2) = [length(z) 1];
        rhoFunction, N2Function, dLnN2Function % function handles
        
        ApU, ApV, ApN
        AmU, AmV, AmN
        A0U, A0V, A0N
        
        UAp, UAm, UA0
        VAp, VAm, VA0
        WAp, WAm
        NAp, NAm, NA0
        
        IMA0, IMAp, IMAm    % InteractionMasks
        EMA0, EMAp, EMAm    % EnergyMasks

        PP, QQ

        t0 = 0
        Ap, Am, A0
        shouldAntiAlias = 1;
        halfK = 0;
        
        offgridModes % subclass should initialize
        ongridModes % This is a cached copy 
        advectionSanityCheck = 0;
        version = 2.1;
    end
    
    properties (Abstract)
        h % all subclasses need to have a function that returns the eigendepths
        
        % These convert the coefficients to their depth integrated energies
%         Apm_HKE_factor
%         Apm_VKE_factor
%         Apm_PE_factor
        Apm_TE_factor
        A0_HKE_factor
        A0_PE_factor
        A0_TE_factor
    end
    
    methods (Abstract)
       u_bar = TransformFromSpatialDomainWithF(self, u)
       w_bar = TransformFromSpatialDomainWithG(self, w)
       u = TransformToSpatialDomainWithF(self, u_bar)
       w = TransformToSpatialDomainWithG(self, w_bar )
       [u,ux,uy,uz] = TransformToSpatialDomainWithFAllDerivatives(self, u_bar)
       [w,wx,wy,wz] = TransformToSpatialDomainWithGAllDerivatives(self, w_bar )
       
       % Needed to add and remove internal waves from the model
       ratio = UmaxGNormRatioForWave(self,k0, l0, j0)
    end
    
    properties (Constant)
        g = 9.81;
    end
    
    methods
        %function self = WaveVortexModel(Lxyz, Nxyz, Nklj, rhobar, N2, dLnN2, latitude, rho0)
        function self = WaveVortexModel(dims, n, z, rhobar, N2, dLnN2, nModes, latitude, rho0)
            % rho0 is optional.
            if length(dims) ~=3 || length(n) ~= 3
                error('The dims and n variables must be of length 3. You need to specify x,y,z');
            end
                                    
            self.Lx = dims(1);
            self.Ly = dims(2);
            self.Lz = dims(3);
            
            self.Nx = n(1);
            self.Ny = n(2);
            self.Nz = n(3);
            self.nModes = nModes;

            dx = self.Lx/self.Nx;
            dy = self.Ly/self.Ny;
            
            self.x = dx*(0:self.Nx-1)'; % periodic basis
            self.y = dy*(0:self.Ny-1)'; % periodic basis
            self.z = z;
            
            dk = 1/self.Lx;          % fourier frequency
            self.k = 2*pi*([0:ceil(self.Nx/2)-1 -floor(self.Nx/2):-1]*dk)';
            dl = 1/self.Ly;          % fourier frequency
            self.l = 2*pi*([0:ceil(self.Ny/2)-1 -floor(self.Ny/2):-1]*dl)';
            self.j = (0:(nModes-1))';

            if self.halfK == 1
                self.k( (self.Nx/2+2):end ) = [];
            end

%             if self.shouldAntiAlias == 1
%                 k_alias = 2*max(abs(self.k))/3;
%                 self.k( abs(self.k) >= k_alias) = [];
%                 self.l( abs(self.l) >= k_alias) = []; % yes, k, we want it isotropic
%                 self.j( abs(self.j) >= 2*max(abs(self.j))/3) = [];
%             end

            self.Nk = length(self.k);
            self.Nl = length(self.l);
            self.Nj = length(self.j);

            
            self.rhobar = rhobar;
            self.N2 = N2;
            self.dLnN2 = dLnN2;
            
            self.Nmax = sqrt(max(N2));
            self.f0 = 2 * 7.2921E-5 * sin( latitude*pi/180 );
            self.latitude = latitude;
            if ~exist('rho0','var')
                self.rho0 = 1025;
            else
                self.rho0 = rho0;
            end

            [X,Y,Z] = ndgrid(self.x,self.y,self.z);
            self.X = X; self.Y = Y; self.Z = Z;
            
            % Now set the initial conditions to zero
            self.Ap = zeros(self.Nk,self.Nl,self.nModes);
            self.Am = zeros(self.Nk,self.Nl,self.nModes);
            self.A0 = zeros(self.Nk,self.Nl,self.nModes);  
            
            % Allow all nonlinear interactions
            self.IMA0 = ones(self.Nk,self.Nl,self.nModes);
            self.IMAp = ones(self.Nk,self.Nl,self.nModes);
            self.IMAm = ones(self.Nk,self.Nl,self.nModes);

            % Allow energy fluxes at all modes
            self.EMA0 = ones(self.Nk,self.Nl,self.nModes);
            self.EMAp = ones(self.Nk,self.Nl,self.nModes);
            self.EMAm = ones(self.Nk,self.Nl,self.nModes);

            if self.shouldAntiAlias == 1
                self.disallowNonlinearInteractionsWithAliasedModes();
                self.freezeEnergyOfAliasedModes();
            end
        end

        function wvmX2 = waveVortexModelWithResolution(self,m)
            wvmX2 = WaveVortexModelHydrostatic([self.Lx self.Ly self.Lz],m, self.latitude, self.rhoFunction, 'N2func', self.N2Function, 'dLnN2func', self.dLnN2Function, 'rho0', self.rho0);
            wvmX2.t0 = self.t0;
            if wvmX2.Nx>=self.Nx && wvmX2.Ny >= self.Ny && wvmX2.nModes >= self.nModes
                kIndices = cat(2,1:(self.Nk/2),(wvmX2.Nk-self.Nk/2 + 1):wvmX2.Nk);
                lIndices = cat(2,1:(self.Nl/2),(wvmX2.Nl-self.Nl/2 + 1):wvmX2.Nl);
                wvmX2.Ap(kIndices,lIndices,1:self.nModes) = self.Ap;
                wvmX2.Am(kIndices,lIndices,1:self.nModes) = self.Am;
                wvmX2.A0(kIndices,lIndices,1:self.nModes) = self.A0;
                
                if self.IsAntiAliased == 1
                    fprintf('You appear to be anti-aliased. When increasing the resolution we will shift the anti-alias filter.\n');
                    AntiAliasMask = self.MaskForAliasedModes();
                    wvmX2.IMA0(kIndices,lIndices,1:self.nModes) = self.IMA0 | AntiAliasMask;
                    wvmX2.IMAp(kIndices,lIndices,1:self.nModes) = self.IMAp | AntiAliasMask;
                    wvmX2.IMAm(kIndices,lIndices,1:self.nModes) = self.IMAm | AntiAliasMask;
                    wvmX2.EMA0(kIndices,lIndices,1:self.nModes) = self.EMA0 | AntiAliasMask;
                    wvmX2.EMAp(kIndices,lIndices,1:self.nModes) = self.EMAp | AntiAliasMask;
                    wvmX2.EMAm(kIndices,lIndices,1:self.nModes) = self.EMAm | AntiAliasMask;

                    wvmX2.disallowNonlinearInteractionsWithAliasedModes();
                    wvmX2.freezeEnergyOfAliasedModes();
                else
                    fprintf('You do NOT appear to be anti-aliased. Thus the interaction masks will be copied as-is.\n');
                    wvmX2.IMA0(kIndices,lIndices,1:self.nModes) = self.IMA0;
                    wvmX2.IMAp(kIndices,lIndices,1:self.nModes) = self.IMAp;
                    wvmX2.IMAm(kIndices,lIndices,1:self.nModes) = self.IMAm;
                    wvmX2.EMA0(kIndices,lIndices,1:self.nModes) = self.EMA0;
                    wvmX2.EMAp(kIndices,lIndices,1:self.nModes) = self.EMAp;
                    wvmX2.EMAm(kIndices,lIndices,1:self.nModes) = self.EMAm;
                end
            else
                error('Reducing resolution not yet implemented. Go for it though, it should be easy.');
            end
        end
        
        function Kh = Kh(self)
            [K,L,~] = ndgrid(self.k,self.l,self.j);
            Kh = sqrt(K.*K + L.*L);
        end 
        
        function Omega = Omega(self)
            [K,L,~] = ndgrid(self.k,self.l,self.j);
            Omega = sqrt(self.g*self.h.*(K.*K + L.*L) + self.f0*self.f0);
        end
        

        
        function rebuildTransformationMatrices(self)
            self.BuildTransformationMatrices(self.PP,self.QQ);
        end

        function self = BuildTransformationMatrices(self,PP,QQ)
            % Build wavenumbers
            [K,L,J] = ndgrid(self.k,self.l,self.j);
            alpha = atan2(L,K);
            K2 = K.*K + L.*L;
            Kh = sqrt(K2);      % Total horizontal wavenumber
            
            f = self.f0;
            g_ = 9.81;
            
            omega = self.Omega;
            if abs(self.f0) < 1e-14 % This handles the f=0 case.
                omega(omega == 0) = 1;
            end
            fOmega = f./omega;
            
            self.PP = PP;
            self.QQ = QQ;
            MakeHermitian = @(f) WaveVortexModel.MakeHermitian(f);
            
            self.iOmega = MakeHermitian(sqrt(-1)*omega);

            if ~exist("PP","var") || isempty(PP)
                PP = ones(size(K));
            end
            if ~exist("QQ","var") || isempty(QQ)
                QQ = ones(size(K));
            end

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
            self.A0U = sqrt(-1)*self.h.*(fOmega./omega) .* L;
            self.A0V = -sqrt(-1)*self.h.*(fOmega./omega) .* K;
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
            
            % Now make the Hermitian conjugate match.
            self.ApU = (1./PP) .* MakeHermitian(self.ApU);
            self.ApV = (1./PP) .* MakeHermitian(self.ApV);
            self.ApN = (1./QQ) .* MakeHermitian(self.ApN);
           
            self.AmU = (1./PP) .* MakeHermitian(self.AmU);
            self.AmV = (1./PP) .* MakeHermitian(self.AmV);
            self.AmN = (1./QQ) .* MakeHermitian(self.AmN);
           
            self.A0U = (1./PP) .* MakeHermitian(self.A0U);
            self.A0V = (1./PP) .* MakeHermitian(self.A0V);
            self.A0N = (1./QQ) .* MakeHermitian(self.A0N);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Transform matrices (Ap,Am,A0) -> (U,V,W,N)
            % These can be pulled from equation C4 in the manuscript
            self.UAp = (cos(alpha)-sqrt(-1)*fOmega.*sin(alpha));
            self.UAm = (cos(alpha)+sqrt(-1)*fOmega.*sin(alpha));
            self.UA0 = -sqrt(-1)*(g_/f)*L;
            
            self.VAp = (sin(alpha)+sqrt(-1)*fOmega.*cos(alpha));
            self.VAm = (sin(alpha)-sqrt(-1)*fOmega.*cos(alpha));
            self.VA0 = sqrt(-1)*(g_/f)*K;
                
            self.WAp = -sqrt(-1)*Kh.*self.h;
            self.WAm = -sqrt(-1)*Kh.*self.h;
            
            self.NAp = -Kh.*self.h./omega;
            self.NAm = Kh.*self.h./omega;
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
            self.UAp = PP .* MakeHermitian(self.UAp);
            self.UAm = PP .* MakeHermitian(self.UAm);
            self.UA0 = PP .* MakeHermitian(self.UA0);
            
            self.VAp = PP .* MakeHermitian(self.VAp);
            self.VAm = PP .* MakeHermitian(self.VAm);
            self.VA0 = PP .* MakeHermitian(self.VA0);
            
            self.WAp = QQ .* MakeHermitian(self.WAp);
            self.WAm = QQ .* MakeHermitian(self.WAm);
            
            self.NAp = QQ .* MakeHermitian(self.NAp);
            self.NAm = QQ .* MakeHermitian(self.NAm);
            self.NA0 = QQ .* MakeHermitian(self.NA0);
        end
          
        function [Ap,Am,A0] = TransformUVEtaToWaveVortex(self,U,V,N,t)
            % This is the 'S^{-1}' operator (C5) in the manuscript
            Ubar = self.TransformFromSpatialDomainWithF(U);
            Vbar = self.TransformFromSpatialDomainWithF(V);
            Nbar = self.TransformFromSpatialDomainWithG(N);
            
            Ap = self.ApU.*Ubar + self.ApV.*Vbar + self.ApN.*Nbar;
            Am = self.AmU.*Ubar + self.AmV.*Vbar + self.AmN.*Nbar;
            A0 = self.A0U.*Ubar + self.A0V.*Vbar + self.A0N.*Nbar;

            if nargin == 5
                phase = exp(-self.iOmega*(t-self.t0));
                Ap = Ap .* phase;
                Am = Am .* conj(phase);
            end
        end
        
        function [U,V,W,N] = TransformWaveVortexToUVWEta(self,Ap,Am,A0,t)
            if nargin == 5
                phase = exp(self.iOmega*(t-self.t0));
                Ap = Ap .* phase;
                Am = Am .* conj(phase);
            end

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
            nNL = u.*DiffFourier(self.x,eta,1,1) + v.*DiffFourier(self.y,eta,1,2) + w.*(DiffSine(self.z,eta,1,3) + eta.*self.dLnN2);
        end
        
        F = NonlinearFluxWithParticlesAtTimeArray(self,t,Y0);
        [Fp,Fm,F0,u,v,w] = NonlinearFluxWithParticlesAtTime(self,t,Ap,Am,A0,x,y,z);
        
        F = NonlinearFluxWithFloatsAndDriftersAtTimeArray(self,t,Y0,z_d);
        [Fp,Fm,F0,u_f,v_f,w_f,u_d,v_d] = NonlinearFluxWithFloatsAndDriftersAtTime(self,t,Ap,Am,A0,x_f,y_f,z_f,x_d,y_d,z_d);
        
        function F = NonlinearFluxAtTimeArray(self,t,Y0)
            [Fp,Fm,F0] = self.NonlinearFluxAtTime(t,Y0{1},Y0{2},Y0{3});
            F = {Fp,Fm,F0};
        end
        
        function [Fp,Fm,F0] = NonlinearFluxAtTime(self,t,Ap,Am,A0)
            % Apply operator T_\omega---defined in (C2) in the manuscript

            % We also apply the interaction masks (IMA*) and energy masks
            % (EMA*). These *could* be precomputed multiplied directly into
            % the coefficients. A quick speed test shows this about a 5%
            % performance hit by not doing this.
            phase = exp(self.iOmega*(t-self.t0));
            Ap = self.IMAp .* Ap .* phase;
            Am = self.IMAm .* Am .* conj(phase);
            A0 = self.IMA0 .* A0;

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
            nNL = -U.*ETAx - V.*ETAy - W.*(ETAz + ETA.*shiftdim(self.dLnN2,-2));
            
            % Now apply the operator S^{-1} and then T_\omega^{-1}
            uNLbar = self.TransformFromSpatialDomainWithF(uNL);
            vNLbar = self.TransformFromSpatialDomainWithF(vNL);
            nNLbar = self.TransformFromSpatialDomainWithG(nNL);

            Fp = self.EMAp .* (self.ApU.*uNLbar + self.ApV.*vNLbar + self.ApN.*nNLbar) .* conj(phase);
            Fm = self.EMAm .* (self.AmU.*uNLbar + self.AmV.*vNLbar + self.AmN.*nNLbar) .* phase;
            F0 = self.EMA0 .* (self.A0U.*uNLbar + self.A0V.*vNLbar + self.A0N.*nNLbar);
        end

        function [Ep,Em,E0] = EnergyFluxAtTime(self,t,Ap,Am,A0)
            [Fp,Fm,F0] = self.NonlinearFluxAtTime(t,Ap,Am,A0);
            % The phase is tricky here. It is wound forward for the flux,
            % as it should be... but then it is wound back to zero. This is
            % equivalent ignoring the phase below here.
            Ep = 2*self.Apm_TE_factor.*real( Fp .* conj(Ap) );
            Em = 2*self.Apm_TE_factor.*real( Fm .* conj(Am) );
            E0 = 2*self.A0_TE_factor.*real( F0 .* conj(A0) );
        end
        
        function [Ep,Em,E0] = EnergyFluxAtTimeInitial(self,t,deltaT,Ap,Am,A0)
      	    [Fp,Fm,F0] = self.NonlinearFluxAtTime(t,Ap,Am,A0);
      	    % The phase is tricky here. It is wound forward for the flux,
      	    % as it should be... but then it is wound back to zero. This is
            % equivalent ignoring the phase below here.

    	    % This equation is C17 in the manuscript, but with addition of 1st term
    	    % on LHS of C16 converted to energy using Apm_TE_factor or A0_TE_factor

    	    % This differs from EnergyFluxAtTime due to the importance of the
    	    % 2*F*F*deltaT in equation C16 at the initial condition.

    	    Ep = 2*Fp.*conj(Fp).*self.Apm_TE_factor*deltaT + 2*self.Apm_TE_factor.*real( Fp .* conj(Ap) );
      	    Em = 2*Fm.*conj(Fm).*self.Apm_TE_factor*deltaT + 2*self.Apm_TE_factor.*real( Fm .* conj(Am) );
      	    E0 = 2*F0.*conj(F0).*self.A0_TE_factor*deltaT + 2*self.A0_TE_factor.*real( F0 .* conj(A0) );
    	end

        [Fp,Fm,F0] = NonlinearFluxForFlowConstituentsAtTime(self,t,Ap,Am,A0,Uconstituent,gradUconstituent)
        [Ep,Em,E0] = EnergyFluxForFlowConstituentsAtTime(self,t,Ap,Am,A0,Uconstituent,gradUconstituent);

    	function [Ep,Em,E0] = EnergyFluxForFlowConstituentsAtTimeInitial(self,t,deltaT,Ap,Am,A0,Uconstituent,gradUconstituent)
    	    [Fp,Fm,F0] = self.NonlinearFluxForFlowConstituentsAtTime(t,Ap,Am,A0,Uconstituent,gradUconstituent);
    	    % The phase is tricky here. It is wound forward for the flux,
    	    % as it should be... but then it is wound back to zero. This is
            % equivalent ignoring the phase below here.

            % This equation is C17 in the manuscript, but with addition of 1st term
            % on LHS of C16 converted to energy using Apm_TE_factor or A0_TE_factor

            % This differs from EnergyFluxForFlowConstituentsAtTime due to the importance of the
            % 2*F*F*deltaT in equation C16 at the initial condition.

            Ep = 2*Fp.*conj(Fp).*self.Apm_TE_factor*deltaT + 2*self.Apm_TE_factor.*real( Fp .* conj(Ap) );
            Em = 2*Fm.*conj(Fm).*self.Apm_TE_factor*deltaT + 2*self.Apm_TE_factor.*real( Fm .* conj(Am) );
            E0 = 2*F0.*conj(F0).*self.A0_TE_factor*deltaT + 2*self.A0_TE_factor.*real( F0 .* conj(A0) );
    	end
 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Reduced interaction models
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function allowNonlinearInteractionsWithModes(self,Ap,Am,A0)
            self.IMA0 = or(self.IMA0,A0);
            self.IMAm = or(self.IMAm,Am);
            self.IMAp = or(self.IMAp,Ap);
        end

        function allowNonlinearInteractionsWithConstituents(self,constituents)
            [ApmMask,A0Mask] = self.MasksForFlowContinuents(constituents);
            self.IMA0 = self.IMA0 | A0Mask;
            self.IMAm = self.IMAm | ApmMask;
            self.IMAp = self.IMAp | ApmMask;
        end

        function disallowNonlinearInteractionsWithConstituents(self,constituents)
            [ApmMask,A0Mask] = self.MasksForFlowContinuents(constituents);
            self.IMA0 = self.IMA0 & ~A0Mask;
            self.IMAm = self.IMAm & ~ApmMask;
            self.IMAp = self.IMAp & ~ApmMask;
        end

        function disallowNonlinearInteractionsWithAliasedModes(self)
            % Uses the 2/3 rule to prevent aliasing of Fourier modes.
            % The reality is that the vertical modes will still alias.
            % http://helper.ipam.ucla.edu/publications/mtws1/mtws1_12187.pdf
            self.shouldAntiAlias = 1;
            
            AntiAliasMask = self.MaskForAliasedModes();
            self.IMA0 = self.IMA0 & ~AntiAliasMask;
            self.IMAm = self.IMAm & ~AntiAliasMask;
            self.IMAp = self.IMAp & ~AntiAliasMask;
        end

        function unfreezeEnergyOfConstituents(self,constituents)
            [ApmMask,A0Mask] = self.MasksForFlowContinuents(constituents);
            self.EMA0 = self.EMA0 | A0Mask;
            self.EMAm = self.EMAm | ApmMask;
            self.EMAp = self.EMAp | ApmMask;
        end

        function freezeEnergyOfConstituents(self,constituents)
            [ApmMask,A0Mask] = self.MasksForFlowContinuents(constituents);
            self.EMA0 = self.EMA0 & ~A0Mask;
            self.EMAm = self.EMAm & ~ApmMask;
            self.EMAp = self.EMAp & ~ApmMask;
        end
        
        function freezeEnergyOfAliasedModes(self)
            % In addition to disallowing interaction to occur between modes
            % that are aliased, you may actually want to disallow energy to
            % even enter the aliased modes.
            AntiAliasMask = self.MaskForAliasedModes();
            self.EMA0 = self.EMA0 & ~AntiAliasMask;
            self.EMAm = self.EMAm & ~AntiAliasMask;
            self.EMAp = self.EMAp & ~AntiAliasMask;
        end

        function clearEnergyFromAliasedModes(self)
            % In addition to disallowing interaction to occur between modes
            % that are aliased, you may actually want to disallow energy to
            % even enter the aliased modes.
            AntiAliasMask = self.MaskForAliasedModes();
            self.A0 = self.A0 .* ~AntiAliasMask;
            self.Am = self.Am .* ~AntiAliasMask;
            self.Ap = self.Ap .* ~AntiAliasMask;
        end

        function flag = IsAntiAliased(self)
            AntiAliasMask = self.MaskForAliasedModes();

            % check if there are zeros at all the anti-alias indices
            flag = all((~self.IMA0 & AntiAliasMask) == AntiAliasMask,'all');
            flag = flag & all((~self.IMAm & AntiAliasMask) == AntiAliasMask,'all');
            flag = flag & all((~self.IMAp & AntiAliasMask) == AntiAliasMask,'all');
            flag = flag & all((~self.EMA0 & AntiAliasMask) == AntiAliasMask,'all');
            flag = flag & all((~self.EMAp & AntiAliasMask) == AntiAliasMask,'all');
            flag = flag & all((~self.EMAm & AntiAliasMask) == AntiAliasMask,'all');
        end

        
        function [omega,k,l] = addForcingWaveModes(self,kModes,lModes,jModes,phi,U,signs)
            [omega,k,l,kIndex,lIndex,jIndex,signIndex] = self.AddGriddedWavesWithWavemodes(kModes,lModes,jModes,phi,U,signs);
            for iMode=1:length(kModes)
                if (signIndex(iMode) == 1 || (kIndex(iMode) == 1 && lIndex(iMode) == 1) )
                    self.EMAp( kIndex(iMode),lIndex(iMode),jIndex(iMode)) = 0;
                    self.EMAp = WaveVortexModel.MakeHermitian(self.EMAp);
                end

                if (signIndex(iMode) == -1 || (kIndex(iMode) == 1 && lIndex(iMode) == 1) )
                    self.EMAm( kIndex(iMode),lIndex(iMode),jIndex(iMode)) = 0;
                    self.EMAm = WaveVortexModel.MakeHermitian(self.EMAm);
                end
            end
        end

        function stirWithConstituents(self,constituents)

        end



        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Energetics (total)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function energy = totalEnergy(self)
            [u,v,w,eta] = self.VariableFieldsAtTime(0,'u','v','w','eta');
            energy = trapz(self.z,mean(mean( u.^2 + v.^2 + w.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
        end
        
        function energy = totalSpectralEnergy(self)
            %             energy = self.inertialEnergy + self.waveEnergy + self.geostrophicEnergy;
            App = self.Ap; Amm = self.Am; A00 = self.A0;
            energy = sum(sum(sum( self.Apm_TE_factor.*( App.*conj(App) + Amm.*conj(Amm) ) + self.A0_TE_factor.*( A00.*conj(A00) ) )));
        end
        
        function energy = totalHydrostaticEnergy(self)
            [u,v,eta] = self.VariableFieldsAtTime(0,'u','v','eta');
            energy = trapz(self.z,mean(mean( u.^2 + v.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
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
                
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Add and remove internal waves from the model
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        period = InitializeWithPlaneWave(self, k0, l0, j0, UAmp, sign)
        
        RemoveAllGriddedWaves(self)
        
        [omega,k,l] = SetGriddedWavesWithWavemodes(self, kMode, lMode, jMode, phi, Amp, signs)
        
        [omega,k,l,kIndex,lIndex,jIndex,signIndex] = AddGriddedWavesWithWavemodes(self, kMode, lMode, jMode, phi, Amp, signs)  
        
        [omega, alpha, k, l, mode, phi, A, norm] = WaveCoefficientsFromGriddedWaves(self);
        
        InitializeWithGMSpectrum(self, GMAmplitude, varargin);
        
        [GM3Dint,GM3Dext] = InitializeWithSpectralFunction(self, GM2D_int, varargin)   ;
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Add and remove geostrophic features from the model
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        SetGeostrophicStreamfunction(self,psi);

        SetInertialMotions(self,u,v);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Eulerian---the dynamical fields on the grid at a given time
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Primary method for accessing the dynamical variables on the
        % internal grid. Valid variable options are 'u', 'v', 'w',
        % 'rho_prime', and 'eta'.
        [varargout] = VariableFieldsAtTime(self, t, varargin);
        
        % Return the velocity field, which is the sum of the gridded and
        % external/free waves at time t. Note that if you do not need w,
        % don't request it and it won't be computed.
        [u,v,w] = VelocityFieldAtTime(self, t);
        
        % Return the density field, which is the sum of the density
        % mean field (variable in z) and the perturbation field
        % (variable in time and space).
        function rho = DensityFieldAtTime(self, t)
            rho = reshape(self.rhobar,1,1,[]) + self.VariableFieldsAtTime(t,'rho_prime');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Lagrangian---return the dynamical fields at a given location and time
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Primary method for accessing the dynamical variables on the at
        % any position or time.
        %
        % The method argument specifies how off-grid values should be
        % interpolated. Use 'exact' for the slow, but accurate, spectral
        % interpolation. Otherwise use 'spline' or some other method used
        % by Matlab's interp function.
        %
        % Valid variable options are 'u', 'v', 'w', 'rho_prime', and
        % 'eta'.
        [varargout] = VariablesAtTimePosition(self,t,x,y,z,interpolationMethod,varargin);
        
        % useful for integration methods where dy/dt is best given with
        % y as a single variable.
        function [u] = VelocityAtTimePositionVector(self,t,p, interpolationMethod,shouldUseW)
            % Return the velocity at time t and position p, where size(p) =
            % [3 n]. The rows of p represent [x,y,z].
            % Optional argument is passed to the interpolation method. I
            % recommend either linear or spline.
            psize = size(p);
            if psize(1) == 3
                [u,v,w] = self.VelocityAtTimePosition(t,p(1,:),p(2,:),p(3,:), interpolationMethod);
                if exist('shouldUseW','var')
                    u = cat(1,u,v,shouldUseW.*w);
                else
                    u = cat(1,u,v,w);
                end
            else
                [u,v,w] = self.VelocityAtTimePosition(t,p(:,1),p(:,2),p(:,3), interpolationMethod);
                if exist('shouldUseW','var')
                    u = cat(2,u,v,shouldUseW.*w);
                else
                    u = cat(2,u,v,w);
                end
            end
        end
        
        function [u] = DrifterVelocityAtTimePositionVector(self,t,p, interpolationMethod)
            % Return the velocity at time t and position p, where size(p) =
            % [3 n]. The rows of p represent [x,y,z].
            psize = size(p);
            if psize(1) == 3
                [u,v,~] = self.VelocityAtTimePosition(t,p(1,:),p(2,:),p(3,:), interpolationMethod);
                u = cat(1,u,v,zeros(1,psize(2)));
            else
                [u,v,~] = self.VelocityAtTimePosition(t,p(:,1),p(:,2),p(:,3), interpolationMethod);
                u = cat(2,u,v,zeros(psize(1),1));
            end
        end
        
        function [u,v,w] = VelocityAtTimePosition(self,t,x,y,z,interpolationMethod)
            if nargout == 3
                [u,v,w] = self.VariablesAtTimePosition(t,x,y,z,interpolationMethod,'u','v','w');
            else
                [u,v] = self.VariablesAtTimePosition(t,x,y,z,interpolationMethod,'u','v');
            end
        end
        
        function rho = DensityAtTimePosition(self,t,x,y,z,interpolationMethod)
            rho = self.rhoFunction(z) + self.VariablesAtTimePosition(t,x,y,z,interpolationMethod,'rho_prime');
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Add and remove off-grid internal waves
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        FillOutWaveSpectrum(self,maxTimeGap)
        
        function RemoveAllExternalWaves(self)
            self.offgridModes.RemoveAllExternalWaves();
        end
        
        function omega = SetExternalWavesWithWavenumbers(self, k, l, j, phi, A, norm)
            omega = self.offgridModes.SetExternalWavesWithWavenumbers(k, l, j, phi, A, norm);
        end
        
        function omega = AddExternalWavesWithWavenumbers(self, k, l, j, phi, A, norm)
            omega = self.offgridModes.AddExternalWavesWithWavenumbers(k, l, j, phi, A, norm);
        end
        
        function k = SetExternalWavesWithFrequencies(self, omega, alpha, j, phi, A, norm)
            k = self.offgridModes.SetExternalWavesWithFrequencies(omega, alpha, j, phi, A, norm);
        end
        
        function k = AddExternalWavesWithFrequencies(self, omega, alpha, j, phi, A, norm)
            k = self.offgridModes.AddExternalWavesWithFrequencies(omega, alpha, j, phi, A, norm);
        end
        
        function [varargout] = ExternalVariableFieldsAtTime(self,t,varargin)
            % Returns the external wave modes at the grid points.
            varargout = cell(size(varargin));
            [X,Y,Z] = ndgrid(self.x,self.y,self.z);
            [varargout{:}] = self.offgridModes.ExternalVariablesAtTimePosition(t,reshape(X,[],1),reshape(Y,[],1), reshape(Z,[],1), varargin{:});
            for iArg=1:length(varargout)
                varargout{iArg} = reshape(varargout{iArg},self.Nx,self.Ny,self.Nz);
            end
        end
        
        function [varargout] = ExternalVariablesAtTimePosition(self,t,x,y,z,varargin)
            varargout = cell(size(varargin));
            [varargout{:}] = self.offgridModes.ExternalVariablesAtTimePosition(t,x,y,z, varargin{:});
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Simulations
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        [deltaT,advectiveDT,oscillatoryDT] = TimeStepForCFL(self, cfl, outputInterval);

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Validation and internal unit testing
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % This is S*S^{-1} and therefore returns the values in
        % wave-vortex space. So, C11 represents Ap and should be 1s
        % where we expected Ap solutions to exist.
        [C11,C21,C31,C12,C22,C32,C13,C23,C33] = ValidateTransformationMatrices(self)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Generate a complete set of wave-vortex coefficients with variance at all
        % physically realizable solution states.
        [ApIO,AmIO,ApIGW,AmIGW,A0G,A0G0,A0rhobar] = GenerateRandomFlowState(self)  

        
        [ApmMask,A0Mask] = MasksForFlowContinuents(self,flowConstituents);
        [IO,SGW,IGW,MDA,SG,IG] = MasksForAllFlowConstituents(self);
        AntiAliasMask= MaskForAliasedModes(self);

        [Qk,Ql,Qj] = ExponentialFilter(self,nDampedModes);
    end
    
    methods (Static)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Check if the matrix is Hermitian. Report errors.
        A = CheckHermitian(A)
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Forces a 3D matrix to be Hermitian, except at k=l=0
        A = MakeHermitian(A)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Returns a matrix the same size as A with 1s at the 'redundant'
        % hermiation indices.
        A = RedundantHermitianCoefficients(A)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Returns a matrix the same size as A with 1s at the Nyquist
        % frequencies.
        A = NyquistWavenumbers(A)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Generate a 3D matrix to be Hermitian, except at k=l=0
        A = GenerateHermitianRandomMatrix( size, shouldExcludeNyquist )
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Takes a Hermitian matrix and resets it back to the real amplitude.
        [A,phi,linearIndex] = ExtractNonzeroWaveProperties(Matrix)
    end
        
        
end 



