classdef WVStratificationConstant < WVStratification
    properties (GetAccess=public, SetAccess=protected)
        N0
        h_0

        F_g,G_g
        DCT, iDCT, DST, iDST
    end

    properties (Dependent)
        FinvMatrix
        GinvMatrix
        FMatrix
        GMatrix
    end

    methods
        function self = WVStratificationConstant(Lz, Nz, options)
            % create matrices for hydrostatics in variable stratification
            %
            % To initialize:
            % 1) Pass (Lz,Nz,N0) as a minimum set, and defaults will be
            % chosen.
            % 2) Pass (Lz,Nz,N0,Nj,latitude,rho0) to fully, uniquely
            % specify the transforms
            %
            % - Topic: Initialization
            % - Declaration: wvt = WVTransformHydrostatic(Lxyz, Nxyz, options)
            % - Parameter Lxyz: length of the domain (in meters) in the three coordinate directions, e.g. [Lx Ly Lz]
            % - Parameter Nxyz: number of grid points in the three coordinate directions, e.g. [Nx Ny Nz]
            % - Parameter rho:  (optional) function_handle specifying the density as a function of depth on the domain [-Lz 0]
            % - Parameter stratification:  (optional) function_handle specifying the stratification as a function of depth on the domain [-Lz 0]
            % - Parameter latitude: (optional) latitude of the domain (default is 33 degrees north)
            % - Parameter rho0: (optional) density at the surface z=0 (default is 1025 kg/m^3)
            % - Returns wvt: a new WVTransformHydrostatic instance
            arguments
                Lz (1,1) double {mustBePositive}
                Nz (1,1) double {mustBePositive}
                options.N0 (1,1) double {mustBePositive} = 5.2e-3
                options.Nj (1,1) double {mustBePositive}
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.planetaryRadius (1,1) double = 6.371e6
                options.rotationRate (1,1) double = 7.2921E-5
                options.latitude (1,1) double = 33
                options.g (1,1) double = 9.81
            end
            dz = Lz/(options.Nz-1);
            options.z = dz*(0:(options.Nz-1))' - Lz; % Cosine basis for DCT-I and DST-I
            options.j = (0:(Nj-1))';
            options.rhoFunction = @(z) -(options.N0*options.N0*options.rho0/options.g)*options.z + options.rho0;
            options.N2Function = @(z) options.N0*options.N0*ones(size(options.z));
            superclassOptions = namedargs2cell(options);
            self@WVStratification(Lz,Nz,superclassOptions{:});

            self.z_int = dz*ones(Nz,1);
            self.z_int(1) = dz/2; self.z_int(end) = dz/2;            
            self.N0 = options.N0;
            self.rho0 = options.rho0;

            self.buildVerticalModeProjectionOperators();
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformation matrices
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function Finv = get.FinvMatrix(wvt)
            % transformation matrix $$F^{-1}$$
            %
            % A matrix that transforms a vector from vertical mode space to physical
            % space.
            %
            % - Topic: Operations — Transformations
            % - Declaration: Finv = FinvMatrix(wvt)
            % - Returns Finv: A matrix with dimensions [Nz Nj]
            arguments
                wvt         WVTransform
            end

            Finv = shiftdim(wvt.F_g(:,1),1) .* wvt.iDCT;

        end

        function F = get.FMatrix(wvt)
            % transformation matrix $$F$$
            %
            % A matrix that transforms a vector from physical
            % space to vertical mode space.
            %
            % - Topic: Operations — Transformations
            % - Declaration: F = FMatrix(wvt)
            % - Returns Finv: A matrix with dimensions [Nz Nj]
            arguments
                wvt         WVTransform
            end

            F = wvt.DCT ./ wvt.F_g(:,1);

        end
        function Ginv = get.GinvMatrix(wvt)
            % transformation matrix $$G^{-1}$$
            %
            % A matrix that transforms a vector from vertical mode space to physical
            % space.
            %
            % - Topic: Operations — Transformations
            % - Declaration: Ginv = GinvMatrix(wvt)
            % - Returns Finv: A matrix with dimensions [Nz Nj]
            arguments
                wvt         WVTransform
            end

            Ginv = shiftdim(wvt.G_g(:,1),1) .* wvt.iDST;

        end
        function G = get.GMatrix(wvt)
            % transformation matrix $$G$$
            %
            % A matrix that transforms a vector from physical
            % space to vertical mode space.
            %
            % - Topic: Operations — Transformations
            % - Declaration: G = GMatrix(wvt)
            % - Returns Ginv: A matrix with dimensions [Nz Nj]
            arguments
                wvt         WVTransform
            end

            G = wvt.DST ./ wvt.G_g(:,1);

        end

        function vm = verticalModes(self)
            vm = InternalModesConstantStratification(N0=self.N0, rho0=self.rho0, zIn=[self.z(1) 0], zOut=self.z, latitude=self.latitude);
        end

        function self = buildVerticalModeProjectionOperators(self)
            % Build the transformation matrices
            self.DCT = WVTransformConstantStratification.CosineTransformForwardMatrix(self.Nz);
            self.DCT = self.DCT(1:self.Nj,:); % dump the Nyquist mode
            self.iDCT = WVTransformConstantStratification.CosineTransformBackMatrix(self.Nz);
            self.iDCT = self.iDCT(:,1:self.Nj); % dump the Nyquist mode
            self.DST = WVTransformConstantStratification.SineTransformForwardMatrix(self.Nz);
            self.DST = cat(1,zeros(1,self.Nz),self.DST);
            self.DST = self.DST(1:self.Nj,:);
            self.iDST = WVTransformConstantStratification.SineTransformBackMatrix(self.Nz);
            self.iDST = cat(2,zeros(self.Nz,1),self.iDST);
            self.iDST = self.iDST(:,1:self.Nj);

            % We renormalization the transformation matrices to directly
            % incorporate normalization of the modes and the DFT.
            [~,~,J] = self.kljGrid;
            M = J*pi/self.Lz;
            N = self.N0;
            g_ = 9.81;

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Normalization for the vertical modes
            % This comes from equations B12 in the manuscript.
            signNorm = -2*(mod(J,2) == 1)+1; % equivalent to (-1)^j
            self.F_g = signNorm .* ((self.h_0).*M)*sqrt(2*g_/(self.Lz*N*N));
            self.G_g = signNorm .* sqrt(2*g_/(self.Lz*N*N));
            self.F_g(J==0) = 2; % j=0 mode is a factor of 2 too big in DCT-I
            self.G_g(J==0) = 1; % j=0 mode doesn't exist for G
        end

        function h = get.h_0(self)
            M = reshape(self.j,[],1)*pi/self.Lz;
            h = (1/self.g)*(self.N0*self.N0)./(M.*M);
            h(1) = self.Lz;
        end
    end

    methods (Access=protected)

    end

    methods (Static)

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = WVStratificationConstant.propertyAnnotationsForStratification();
        end
        function vars = classRequiredPropertyNames()
            vars = WVStratificationConstant.requiredPropertiesForStratification();
        end

        function requiredPropertyNames = namesOfRequiredPropertiesForStratification()
            requiredPropertyNames = WVStratification.namesOfRequiredPropertiesForStratification();
            requiredPropertyNames = union(requiredPropertyNames,WVStratificationConstant.newRequiredPropertyNames());
        end

        function newRequiredPropertyNames = newRequiredPropertyNames()
            newRequiredPropertyNames = {'N0'};
        end

        function propertyAnnotations = propertyAnnotationsForStratification()
            % return array of property annotations initialized by default
            %
            % This function returns annotations for all properties of the
            % WVStratificationConstant class (as well as its
            % superclass).
            %
            % - Topic: Internal
            % - Declaration: propertyAnnotations = WVStratificationConstant.propertyAnnotationsForStratification()
            % - Returns propertyAnnotations: array of WVPropertyAnnotation instances
            propertyAnnotations = WVStratification.propertyAnnotationsForStratification();
            propertyAnnotations(end+1) = CANumericProperty('N0',{},'rad s^{-1}', 'buoyancy frequency of the no-motion density');
        end

        function [Lz,Nz,options] = requiredPropertiesForStratificationFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                Lz (1,1) double {mustBePositive}
                Nz (1,1) double {mustBePositive}
                options
            end
            [Lz,Nz,stratOptions] = WVStratification.requiredPropertiesForStratificationFromGroup(group);
            vars = CAAnnotatedClass.propertyValuesFromGroup(group,WVStratificationConstant.newRequiredPropertyNames);
            newOptions = namedargs2cell(vars);
            options = cat(2,stratOptions,newOptions);
        end

        function stratification = stratificationFromFile(path)
            arguments (Input)
                path char {mustBeNonempty}
            end
            arguments (Output)
                stratification WVStratificationVariable {mustBeNonempty}
            end
            ncfile = NetCDFFile(path);
            stratification = WVStratificationConstant.stratificationFromGroup(ncfile);
        end

        function stratification = stratificationFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                stratification WVStratificationVariable {mustBeNonempty}
            end
            CAAnnotatedClass.throwErrorIfMissingProperties(group,WVStratificationConstant.namesOfRequiredPropertiesForStratification);
            [Lz,Nz,options] = WVStratificationConstant.requiredPropertiesForStratificationFromGroup(group);
            stratification = WVStratificationConstant(Lz,Nz,options{:});
        end
        
    end
end