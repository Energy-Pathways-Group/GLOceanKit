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
        Fp % size(F)=[Nz x Nmodes+1], barotropic mode AND extra Nyquist mode
        Gp % size(G)=[Nz-2 x Nmodes-1], no barotropic mode
        h % [1 x Nmodes]
        
        Pf % Preconditioner for F, size(Fp)=[1 Nmodes+1]. u = F*m, u = Pf*inv(Pf)*F*m, so Fp==inv(Pf)*F. 
        Pg % Preconditioner for G, size(Gp)=[1 Nmodes-1]. eta = G*m, eta = Pg*inv(Pg)*G*m, so Gp==inv(Pg)*G. 

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
            z = linspace(-dims(3),0,nModes*10)';
            im = InternalModesWKBSpectral(rhoFunc,[-dims(3) 0],z,latitude);
            im.normalization = Normalization.kConstant;
            im.upperBoundary = UpperBoundary.rigidLid;
            z = im.GaussQuadraturePointsForModesAtFrequency(nModes+1,0);
            
            % Now go compute the appropriate number of modes at the
            % quadrature points.
            im = InternalModesWKBSpectral(rhoFunc,[-dims(3) 0],z,latitude,'nModes',nModes);
            im.normalization = Normalization.kConstant;
            im.upperBoundary = UpperBoundary.rigidLid;
            
            % This is enough information to initialize
            self@WaveVortexModel(dims, n, z, im.rho, im.N2, nModes, latitude, im.rho_function(0));
            self.rhoFunction = im.rho_function;
            self.N2Function = im.N2_function;
            self.internalModes = im;
            [F,G,self.h] = self.internalModes.ModesAtFrequency(0);

            % Add the barotropic mode
            F = cat(2,ones(self.Nz,1),F);
            self.h = cat(2,1,self.h);
            self.h = shiftdim(self.h,-1);

            % Compute the precondition matrices (really, diagonals)
            self.Pf = max(abs(F),[],1);
            self.Pg = max(abs(G),[],1);

            % Now create the actual transformation matrices
            self.Fp = F./Pf;
            self.Gp = G./Gf;

            self.BuildTransformationMatrices();
            self.offgridModes = WaveVortexModelOffGrid(im,latitude, self.N2Function);
        end
                                
%         function self = BuildTransformationMatrices(self)
%             BuildTransformationMatrices@WaveVortexModel(self);
%             
%             
% 
% 
%             % We renormalization the transformation matrices to directly
%             % incorporate normalization of the modes and the DFT.          
%             [~,~,J] = ndgrid(self.k,self.l,self.j);
%             M = J*pi/self.Lz;
%             N = self.N0;
%             f = self.f0; 
%             g_ = 9.81;
%        
%             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             % Normalization for the vertical modes
%             % This comes from equations B12 in the manuscript.
%             signNorm = -2*(mod(J,2) == 1)+1; % equivalent to (-1)^j
%             self.F = signNorm .* ((self.h).*M)*sqrt(2*g_/(self.Lz*(N*N-f*f)));
%             self.G = signNorm .* sqrt(2*g_/(self.Lz*(N*N-f*f)));
%             self.F(:,:,1) = 2; % j=0 mode is a factor of 2 too big in DCT-I
%             self.G(:,:,1) = 1; % j=0 mode doesn't exist for G
%   
%             % Now make the Hermitian conjugate match.
%             iFTransformScaling = 2./(self.Nx*self.Ny*self.F);
%             iGTransformScaling = 2./(self.Nx*self.Ny*self.G);
%             self.ApU = iFTransformScaling .* self.ApU;
%             self.ApV = iFTransformScaling .* self.ApV;
%             self.ApN = iGTransformScaling .* self.ApN;
%             
%             self.AmU = iFTransformScaling .* self.AmU;
%             self.AmV = iFTransformScaling .* self.AmV;
%             self.AmN = iGTransformScaling .* self.AmN;
%             
%             self.A0U = iFTransformScaling .* self.A0U;
%             self.A0V = iFTransformScaling .* self.A0V;
%             self.A0N = iGTransformScaling .* self.A0N;
%                         
%             % Now make the Hermitian conjugate match AND pre-multiply the
%             % coefficients for the transformations.
%             FTransformScaling = 0.5*self.Nx*self.Ny*self.F;
%             self.UAp = FTransformScaling .* self.UAp;
%             self.UAm = FTransformScaling .* self.UAm;
%             self.UA0 = FTransformScaling .* self.UA0;
%             
%             self.VAp = FTransformScaling .* self.VAp;
%             self.VAm = FTransformScaling .* self.VAm;
%             self.VA0 = FTransformScaling .* self.VA0;
%             
%             GTransformScaling = 0.5*self.Nx*self.Ny*self.G;
%             self.WAp = GTransformScaling .* self.WAp;
%             self.WAm = GTransformScaling .* self.WAm;
%             
%             self.NAp = GTransformScaling .* self.NAp;
%             self.NAm = GTransformScaling .* self.NAm;
%             self.NA0 = GTransformScaling .* self.NA0;
%         end
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %
%         % Energetics
%         %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         function value = get.Apm_TE_factor(self)
%             value = self.h; % factor of 2 larger than in the manuscript
%             value(:,:,1) = self.Lz;
%         end
%         
%         function value = get.A0_HKE_factor(self)
%             [K,L,J] = ndgrid(self.k,self.l,self.j);
%             K2 = K.*K + L.*L;
%             M = J*pi/self.Lz;
%             
%             % This comes from equation (3.10) in the manuscript, but using
%             % the relation from equation A2b
%             % omega = sqrt(self.g*h.*K2 + self.f0*self.f0);
%             % value = (self.g/(self.f0*self.f0)) * (omega.*omega - self.f0*self.f0) .* (self.N0*self.N0 - omega.*omega) / (2 * (self.N0*self.N0 - self.f0*self.f0) );
%             value = (self.g^3/(self.f0*self.f0)) * K2.*self.h.*self.h.*M.*M / (2 * (self.N0*self.N0 - self.f0*self.f0) ); % factor of 2 larger than in the manuscript
%             value(:,:,1) = (self.g^2/(self.f0*self.f0)) * K2(:,:,1) * self.Lz/2;
%         end
%         function value = get.A0_PE_factor(self)
%             value = self.g*self.N0*self.N0/(self.N0*self.N0-self.f0*self.f0)/2; % factor of 2 larger than in the manuscript
%         end
%         function value = get.A0_TE_factor(self)
%             value = self.A0_HKE_factor + self.A0_PE_factor;
%         end
%           
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %
%         % Wave properties
%         %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         
%         function cg_x = get.cg_x(self)
%             [K,L,J] = ndgrid(self.k,self.l,self.j);
%             K2 = K.*K + L.*L;
%             M = J*pi/self.Lz;
%             Omega = sqrt( (self.N0*self.N0*K2+self.f0*self.f0*M.*M)./(K2+M.*M) );
%             cg_x = (K./Omega) .*M.*M .* (self.N0*self.N0-self.f0*self.f0)./(M.*M+K2).^2;
%             cg_x(isnan(cg_x)) = 0;
%         end
%         
%         function cg_y = get.cg_y(self)
%             [K,L,J] = ndgrid(self.k,self.l,self.j);
%             K2 = K.*K + L.*L;
%             M = J*pi/self.Lz;
%             Omega = sqrt( (self.N0*self.N0*K2+self.f0*self.f0*M.*M)./(K2+M.*M) );
%             cg_y = (L./Omega) .* M.*M .* (self.N0*self.N0-self.f0*self.f0)./(M.*M+K2).^2;
%             cg_y(isnan(cg_y)) = 0;
%         end
%         
%         function cg_z = get.cg_z(self)
%             [K,L,J] = ndgrid(self.k,self.l,self.j);
%             K2 = K.*K + L.*L;
%             M = J*pi/self.Lz;
%             Omega = sqrt( (self.N0*self.N0*K2+self.f0*self.f0*M.*M)./(K2+M.*M) );
%             cg_z = -(M./Omega) .* K2 .* (self.N0*self.N0-self.f0*self.f0)./(M.*M+K2).^2;
%             cg_z(isnan(cg_z)) = 0;
%         end
%         
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %
%         % Transformations to and from the spatial domain
%         %
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%       
        function u_bar = TransformFromSpatialDomainWithF(self, u)

        end
        
        function w_bar = TransformFromSpatialDomainWithG(self, w)

        end
        
        function u = TransformToSpatialDomainWithF(self, u_bar)
        end  
                
        function w = TransformToSpatialDomainWithG(self, w_bar )
        end
        
        function [u,ux,uy,uz] = TransformToSpatialDomainWithFAllDerivatives(self, u_bar)
        end  
        
        function [w,wx,wy,wz] = TransformToSpatialDomainWithGAllDerivatives(self, w_bar )
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



