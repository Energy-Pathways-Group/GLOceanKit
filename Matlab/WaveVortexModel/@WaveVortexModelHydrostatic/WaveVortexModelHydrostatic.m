classdef WaveVortexModelHydrostatic < WaveVortexModel
    % 3D hydrostatic Boussinesq model with arbitrary stratification solved
    % in wave-vortex space
    %
    % Couple of different initialization paths:
    % 1) You want to run this as a prognostic model and therefore want
    %    the chebyshev points automatically found for you
    %       Init([Lx Ly Lz], [Nx Ny nModes], latitude, rho)
    %
    % 2) You want to run this as a diagnostic model and therefore want
    %    to specify the depths and modes yourself
    %       Init([Lx Ly Lz], [Nx Ny nModes], latitude, rho, 'zgrid', z)

    properties
        internalModes
        
        % Transformation matrices
        PFinv, QGinv % size(PFinv,PGinv)=[Nz x Nmodes]
        PF, QG % size(PF,PG)=[Nmodes x Nz]
        h % [1 x Nmodes]
        
        P % Preconditioner for F, size(P)=[1 Nmodes]. F*u = uhat, (PF)*u = P*uhat, so ubar==P*uhat 
        Q % Preconditioner for G, size(Q)=[1 Nmodes]. G*eta = etahat, (QG)*eta = Q*etahat, so etabar==Q*etahat. 

        Apm_TE_factor
        A0_HKE_factor
        A0_PE_factor
        A0_TE_factor
    end
        
    methods
         
        function self = WaveVortexModelHydrostatic(dims, n, latitude, rhoFunc, varargin)
            if length(dims) ~=3 || length(n) ~= 3
                error('The dims and n variables must be of length 3. You need to specify x,y,z');
            end
            
            nargs = length(varargin);
            if mod(nargs,2) ~= 0
                error('Arguments must be given as name/value pairs');
            end
            
            % First thing we do is find the Gauss-quadrature points for
            % this stratification profile.
            nModes = n(3);
            Nz = nModes+1;
            z = linspace(-dims(3),0,Nz*10)';
            im = InternalModesWKBSpectral(rhoFunc,[-dims(3) 0],z,latitude);
            im.normalization = Normalization.kConstant;
            im.upperBoundary = UpperBoundary.rigidLid;
            z = im.GaussQuadraturePointsForModesAtFrequency(nModes+1,0);
                        
            % If the user requests nModes---we will have that many fully
            % resolved modes. It works as follows for rigid lid:
            % - There is one barotropic mode that appears in F
            % - There are nModes-1 *internal modes* for G and F.
            % - We compute the nModes+1 internal mode for F, to make it
            % complete.
            % This is nModes+1 grid points necessary to make this happen.
            % This should make sense because there are nModes-1 internal
            % modes, but the boundaries.
            im = InternalModesWKBSpectral(rhoFunc,[-dims(3) 0],z,latitude,'nModes',nModes);
            im.normalization = Normalization.kConstant;
            im.upperBoundary = UpperBoundary.rigidLid;

            % This is enough information to initialize
            self@WaveVortexModel(dims, [n(1) n(2) Nz], z, im.rho, im.N2, nModes, latitude, im.rho_function(0));
            self.rhoFunction = rhoFunc;
            self.N2Function = im.N2_function;
            self.internalModes = im;

            self.BuildTransformationMatrices();
            self.offgridModes = WaveVortexModelOffGrid(im,latitude, self.N2Function);
        end
                                
        function self = BuildTransformationMatrices(self)
            % Now go compute the appropriate number of modes at the
            % quadrature points.
            [Finv,Ginv,self.h] = self.internalModes.ModesAtFrequency(0);
            
            % Make these matrices invertible by adding the barotropic mode
            % to F, and removing the boundaries of G.
            Finv = cat(2,ones(self.Nz,1),Finv);
            Ginv = Ginv(2:end-1,1:end-1);

            % Compute the precondition matrices (really, diagonals)
            self.P = max(abs(Finv),[],1);
            self.Q = max(abs(Ginv),[],1);

            % Now create the actual transformation matrices
            self.PFinv = Finv./self.P;
            self.QGinv = Ginv./self.Q;
            self.PF = inv(self.PFinv);
            self.QG = inv(self.QGinv);

            % size(F)=[Nz x Nmodes+1], barotropic mode AND extra Nyquist mode
            % but, we will only multiply by vectors [Nmodes 1], so dump the
            % last column. Now size(Fp) = [Nz x Nmodes].
            self.PFinv = self.PFinv(:,1:end-1);

            % size(Finv)=[Nmodes+1, Nz], but we don't care about the last mode
            self.PF = self.PF(1:end-1,:);
            
            % size(G) = [Nz-2, Nmodes-1], need zeros for the boundaries
            % and add the 0 barotropic mode, so size(G) = [Nz, Nmodes],
            self.QGinv = cat(2,zeros(self.Nz,1),cat(1,zeros(1,self.nModes-1),self.QGinv,zeros(1,self.nModes-1)));

            % size(Ginv) = [Nmodes-1, Nz-2], need a zero for the barotropic
            % mode, but also need zeros for the boundary
            self.QG = cat(2,zeros(self.nModes,1), cat(1,zeros(1,self.Nz-2),self.QG),zeros(self.nModes,1));

            % want size(h)=[1 1 nModes]
            self.h = cat(2,1,self.h(1:end-1)); % remove the extra mode at the end
            self.h = shiftdim(self.h,-1);

            self.P = shiftdim(self.P(1:end-1),-1);
            self.Q = shiftdim(cat(2,1,self.Q),-1);
            
            PP = self.Nx*self.Ny*self.P;
            QQ = self.Nx*self.Ny*self.Q;

            BuildTransformationMatrices@WaveVortexModel(self,PP,QQ);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Energetics
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function value = get.Apm_TE_factor(self)
            value = repmat(self.h,self.Nx,self.Ny); % factor of 2 larger than in the manuscript
            value(:,:,1) = self.Lz;
        end
        
        function value = get.A0_HKE_factor(self)
            [K,L,~] = ndgrid(self.k,self.l,self.j);
            K2 = K.*K + L.*L;

            value = (self.g^2/(self.f0*self.f0)) * K2 .* self.Apm_TE_factor/2;
        end
        function value = get.A0_PE_factor(self)
            value = self.g*ones(self.Nx,self.Ny,self.nModes)/2;
        end
        function value = get.A0_TE_factor(self)
            value = self.A0_HKE_factor + self.A0_PE_factor;
        end
          
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations to and from the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function u_bar = TransformFromSpatialDomainWithF(self, u)
            % hydrostatic modes commute with the DFT
            u = permute(u,[3 1 2]); % keep adjacent in memory
            u = reshape(u,self.Nz,[]);
            u_bar = self.PF*u;
            u_bar = reshape(u_bar,self.nModes,self.Nx,self.Ny);
            u_bar = permute(u_bar,[2 3 1]);
            u_bar = fft(fft(u_bar,self.Nx,1),self.Ny,2);
        end
        
        function w_bar = TransformFromSpatialDomainWithG(self, w)
            % hydrostatic modes commute with the DFT
            w = permute(w,[3 1 2]); % keep adjacent in memory
            w = reshape(w,self.Nz,[]);
            w_bar = self.QG*w;
            w_bar = reshape(w_bar,self.nModes,self.Nx,self.Ny);
            w_bar = permute(w_bar,[2 3 1]);
            w_bar = fft(fft(w_bar,self.Nx,1),self.Ny,2);
        end
        
        function u = TransformToSpatialDomainWithF(self, u_bar)
            % hydrostatic modes commute with the DFT
            u_bar = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');
            u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
            u_bar = reshape(u_bar,self.nModes,[]);
            u = self.PFinv*u_bar;
            u = reshape(u,self.Nz,self.Nx,self.Ny);
            u = permute(u,[2 3 1]);
        end  
                
        function w = TransformToSpatialDomainWithG(self, w_bar )
            % hydrostatic modes commute with the DFT
            w_bar = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');
            w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
            w_bar = reshape(w_bar,self.nModes,[]);
            w = self.QGinv*w_bar;
            w = reshape(w,self.Nz,self.Nx,self.Ny);
            w = permute(w,[2 3 1]);
        end
        
        function [u,ux,uy,uz] = TransformToSpatialDomainWithFAllDerivatives(self, u_bar)
            u_bar = ifft(ifft(u_bar,self.Nx,1),self.Ny,2,'symmetric');

            u_bar = permute(u_bar,[3 1 2]); % keep adjacent in memory
            u_bar = reshape(u_bar,self.nModes,[]);
            u = self.PFinv*u_bar;
            u = reshape(u,self.Nz,self.Nx,self.Ny);
            u = permute(u,[2 3 1]);

            ux = ifft( sqrt(-1)*self.k.*fft(u,self.Nx,1), self.Nx, 1,'symmetric');
            uy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(u,self.Ny,2), self.Ny, 2,'symmetric');

            uz = self.QGinv*( squeeze(self.Q./self.P).*u_bar );
            uz = reshape(uz,self.Nz,self.Nx,self.Ny);
            uz = permute(uz,[2 3 1]);
            uz = (-shiftdim(self.N2,-2)/self.g).*uz;
        end  
        
        function [w,wx,wy,wz] = TransformToSpatialDomainWithGAllDerivatives(self, w_bar )
            w_bar = ifft(ifft(w_bar,self.Nx,1),self.Ny,2,'symmetric');

            w_bar = permute(w_bar,[3 1 2]); % keep adjacent in memory
            w_bar = reshape(w_bar,self.nModes,[]);
            w = self.QGinv*w_bar;
            w = reshape(w,self.Nz,self.Nx,self.Ny);
            w = permute(w,[2 3 1]);

            wx = ifft( sqrt(-1)*self.k.*fft(w,self.Nx,1), self.Nx, 1,'symmetric');
            wy = ifft( sqrt(-1)*shiftdim(self.l,-1).*fft(w,self.Ny,2), self.Ny, 2,'symmetric');
            
            wz = self.PFinv* ( squeeze(self.P./(self.Q .* self.h)) .* w_bar);
            wz = reshape(wz,self.Nz,self.Nx,self.Ny);
            wz = permute(wz,[2 3 1]);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Needed to add and remove internal waves from the model
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function ratio = UmaxGNormRatioForWave(self,k0, l0, j0)

        end   
        
    end
   
        
        
end 



