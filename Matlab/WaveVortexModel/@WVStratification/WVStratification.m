classdef WVStratification < WVRotatingFPlane
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here
    properties (GetAccess=public, SetAccess=protected)
        % eta_true operation needs rhoFunction
        rhoFunction, N2Function%, dLnN2Function = []
        rho0, rho_nm0, N2%, dLnN2
        z, j
        z_int
        verticalModes

        % length of the z-dimension
        %
        % - Topic: Domain attributes — Spatial grid
        Lz
    end

    properties (Dependent)
        Nz, Nj
    end

    methods (Abstract)
        u_z = diffZF(self,u,options);
        w_z = diffZG(self,w,options);
    end

    properties (Abstract)
        FinvMatrix
        GinvMatrix

        FMatrix
        GMatrix
    end

    methods
        function value = get.Nz(self)
            value=length(self.z);
        end

        function value = get.Nj(self)
            value=length(self.j);
        end

        function flag = isDensityInValidRange(self)
            % checks if the density field is a valid adiabatic re-arrangement of the base state
            %
            % This is probably best re-defined as a dynamical variable.
            %
            % - Topic: Stratification — Validation
            % - Declaration: flag = isDensityInValidRange()
            % - Returns flag: a boolean
            flag = ~(any(self.rho_total(:) < min(self.rho_nm0)) | any(self.rho_total(:) > max(self.rho_nm0)));
        end

        function effectiveVerticalGridResolution = effectiveVerticalGridResolution(self)
            %returns the effective vertical grid resolution in meters
            %
            % The effective grid resolution is the highest fully resolved
            % wavelength in the model. This value takes into account
            % anti-aliasing, and is thus appropriate for setting damping
            % operators.
            %
            % - Topic: Stratification
            % - Declaration: flag = effectiveVerticalGridResolution(other)
            % - Returns effectiveVerticalGridResolution: double
            arguments
                self WVStratification
            end
            effectiveVerticalGridResolution = min(diff(self.z));
        end

        function throwErrorIfDensityViolation(self,options)
            % checks if the proposed coefficients are a valid adiabatic re-arrangement of the base state
            %
            % Given some proposed new set of values for A0, Ap, Am, will
            % the fluid state violate our density condition? If yes, then
            % throw an error and tell the user about it.
            arguments
                self WVStratification
                options.Ap double = 0
                options.Am double = 0
                options.A0 double = 0
                options.additionalErrorInfo = sprintf(' ');
            end
            % Trying to be extra careful here, so we include the wave part
            % of the flow because we do not really know how everything will
            % add together in the end.
            if ~isscalar(options.Am) || ~isscalar(options.Ap)
                rho_total = reshape(self.rho_nm0,1,1,[]) + (self.rho0/self.g) * shiftdim(self.N2,-2) .* self.transformToSpatialDomainWithG(A0=self.NA0.*options.A0,Apm=self.NAp.*options.Ap + self.NAm.*options.Am);
            else
                rho_total = reshape(self.rho_nm0,1,1,[]) + (self.rho0/self.g) * shiftdim(self.N2,-2) .* self.transformToSpatialDomainWithG(A0=self.NA0.*options.A0);
            end
            densityViolation = any(rho_total(:) < min(self.rho_nm0)) | any(rho_total(:) > max(self.rho_nm0));
            if densityViolation == 1
                errorString = sprintf('The no-motion density minus rho0 spans from %.3f kg/m^{3} at the surface to %.3f kg/m^{3} at the bottom. Any adiabatic re-arrangement of the fluid requires the density anomaly stay within this range. ',self.rho_nm0(end)-self.rho0,self.rho_nm0(1)-self.rho0);
                minString = sprintf('\tminimum density: %.3f kg/m^{3}\n',min(rho_total(:))-self.rho0);
                maxString = sprintf('\tmaximum density: %.3f kg/m^{3}\n',max(rho_total(:))-self.rho0);
                errorStruct.message = [errorString,options.additionalErrorInfo,minString,maxString];
                errorStruct.identifier = 'WVTransform:DensityBoundsViolation';
                error(errorStruct);
            end
        end

        function cheb_function = chebfunForZArray(self,my_z_vector)
            cheb_function = chebfun( @(z) interp1(self.z,my_z_vector,z,'spline'),[min(self.z) max(self.z)],'splitting','on');
        end
    end

    methods (Access=protected)
        function self = WVStratification(Lz,Nz,options,rotatingOptions)
            % By assumption, the modes j=0:(Nj-1) --- because for now, we
            % require completeness.
            % 
            % The required variables for a unique specification of the
            % stratification are (Lz,Nz,N2,Nj,latitude,rho0). However, (Nj,
            % latitude, rho0) all have defaults.
            %
            % To initialize:
            % 1) Pass (Lz,Nz,N2) as a minimum set, and defaults will be
            % chosen (or rho instead of N2).
            % 2) You can also pass everything, and override the defaults,
            % but the assumption is that z is the quadrature grid.
            % 3) verticalModes is strictly optional, and will be created
            % for you if you do not pass it.
            %
            % At the end of initialization, the verticalModes class is
            % initialized with the quadrature point
            arguments
                Lz (1,1) double {mustBePositive}
                Nz (1,1) double {mustBePositive}
                options.z (:,1) double {mustBeNonempty} % quadrature points!
                options.j (:,1) double {mustBeNonempty}
                options.Nj (1,1) double {mustBePositive}
                options.rhoFunction function_handle = @isempty
                options.N2Function function_handle = @isempty
                options.rho0 (1,1) double {mustBePositive} = 1025
                rotatingOptions.planetaryRadius (1,1) double = 6.371e6
                rotatingOptions.rotationRate (1,1) double = 7.2921E-5
                rotatingOptions.latitude (1,1) double = 33
                rotatingOptions.g (1,1) double = 9.81
            end

            if isequal(options.N2Function,@isempty) && isequal(options.rhoFunction,@isempty)
                error('You must specify either rho or N2.');
            end
            % For rigid lid:
            % - There is one barotropic mode that appears in F
            % - There are nModes-1 *internal modes* for G and F.
            % - We compute the nModes+1 internal mode for F, to make it
            % complete.
            % This is nModes+1 grid points necessary to make this happen.
            % This should make sense because there are nModes-1 internal
            % modes, but the boundaries.
            optionCell = namedargs2cell(rotatingOptions);
            self@WVRotatingFPlane(optionCell{:});

            self.z=options.z;
            self.Lz = Lz;

            if isfield(options,'j')
                self.j=options.j;
            else
                if isfield(options,'Nj')
                    Nj = options.Nj;
                else
                    Nj = self.Nz-1;
                end
                self.j = (0:(Nj-1))';
            end

            self.rho0 = options.rho0;
            self.N2 = options.N2Function(self.z);
            self.N2Function = options.N2Function;
            self.rhoFunction = options.rhoFunction;
            self.rho_nm0 = self.rhoFunction(self.z);
        end

        function [P,Q,PFinv,PF,QGinv,QG,h,w] = verticalProjectionOperatorsForGeostrophicModes(self,Nj)
            % Now go compute the appropriate number of modes at the
            % quadrature points.
            self.verticalModes.normalization = Normalization.geostrophic;
            self.verticalModes.upperBoundary = UpperBoundary.rigidLid;
            [Finv,Ginv,h] = self.verticalModes.ModesAtFrequency(0);
            [P,Q,PFinv,PF,QGinv,QG,h,w] = WVStratification.verticalProjectionOperatorsWithRigidLid(Finv,Ginv,h,Nj,self.verticalModes.Lz);
        end

        function [P,Q,PFinv,PF,QGinv,QG,h] = verticalProjectionOperatorsForIGWModes(self,k,Nj)
            % Now go compute the appropriate number of modes at the
            % quadrature points.
            self.verticalModes.normalization = Normalization.kConstant;
            self.verticalModes.upperBoundary = UpperBoundary.rigidLid;
            [Finv,Ginv,h] = self.verticalModes.ModesAtWavenumber(k);
            [P,Q,PFinv,PF,QGinv,QG,h] = WVStratification.verticalProjectionOperatorsWithRigidLid(Finv,Ginv,h,Nj,self.verticalModes.Lz);
        end
    end

    methods (Static)
        z = quadraturePointsForStratifiedFlow(Lz,Nz,options);
        [P,Q,PFinv,PF,QGinv,QG,h,w] = verticalProjectionOperatorsWithRigidLid(Finv,Ginv,h,Nj,Lz);
        [P,Q,PFinv,PF,QGinv,QG,h] = verticalProjectionOperatorsWithFreeSurface(Finv,Ginv,h,Nj,Lz);
    end
    methods (Static, Hidden=true)
        % All the metadata has to be defined at the class level---so static
        % methods. Thus we have,
        % -dimensionAnnotationsForStratification
        % -propertyAnnotationsForStratification
        % -methodAnnotationsForStratification
        % -requiredDimensionsForStratification
        % -requiredVariablesForStratification
        % function vars = requiredVariablesForStratification()
        %     vars = {'rho0','rho_nm0','N2'};
        % end
        % 
        % function dims = requiredDimensionsForStratification()
        %     dims = {'z','j'};
        % end


        function Nz = verticalResolutionForHorizontalResolution(Lxy,Lz,Nx,options,rotatingOptions)
            % optimal resolution for a QG model, such that k_max =
            % 1/\lambda_max. WKB scaling would suggest
            % Nz = Nxy * (pi * sqrt(Lr2(2)))/Lxy
            arguments
                Lxy (1,1) double {mustBePositive}
                Lz (1,1) double {mustBePositive}
                Nx (1,1) double {mustBePositive}
                options.rhoFunction function_handle = @isempty
                options.N2Function function_handle = @isempty
                options.rho0 (1,1) double {mustBePositive} = 1025
                rotatingOptions.rotationRate (1,1) double = 7.2921E-5
                rotatingOptions.latitude (1,1) double = 33
                rotatingOptions.g (1,1) double = 9.81
            end
            z = linspace(-Lz,0,500);
            nEVP = max(256,Nx);
            nModes = floor(nEVP/2.1);
            im = InternalModesWKBSpectral(N2=options.N2Function,zIn=[-Lz 0],zOut=z,latitude=rotatingOptions.latitude,rho0=options.rho0,nModes=nModes,nEVP=nEVP,rotationRate=rotatingOptions.rotationRate,g=rotatingOptions.g);
            im.normalization = Normalization.geostrophic;
            im.upperBoundary = UpperBoundary.rigidLid;
            [~,~,h] = im.ModesAtFrequency(0);
            Lr = sqrt(rotatingOptions.g*h)/im.f0;
            kmax = pi*Nx/Lxy;

            [~,Nz] = min( abs(Lr*kmax-1) );
        end

        function requiredPropertyNames = namesOfRequiredPropertiesForStratification()
            requiredPropertyNames = WVRotatingFPlane.namesOfRequiredPropertiesForRotatingFPlane();
            requiredPropertyNames = union(requiredPropertyNames,{'z','Lz','j','rho0','N2Function'});
        end

        function propertyAnnotations = propertyAnnotationsForStratification()
            % return array of CAPropertyAnnotations initialized by default
            %
            % This function returns annotations for all properties of the
            % WVStratification class.
            %
            % - Topic: Developer
            % - Declaration: propertyAnnotations = WVStratification.propertyAnnotationsForStratification()
            % - Returns propertyAnnotations: array of CAPropertyAnnotation instances
            propertyAnnotations = WVRotatingFPlane.propertyAnnotationsForRotatingFPlane();

            propertyAnnotations(end+1) = CADimensionProperty('z', 'm', 'z coordinate');
            propertyAnnotations(end).attributes('standard_name') = 'height_above_mean_sea_level';
            propertyAnnotations(end).attributes('positive') = 'up';
            propertyAnnotations(end).attributes('axis') = 'Z';

            propertyAnnotations(end+1) = CADimensionProperty('j', 'mode number', 'vertical mode number');
            propertyAnnotations(end+1) = CANumericProperty('Nj',{},'', 'points in the j-coordinate, `length(z)`', detailedDescription='- topic: Domain Attributes — Grid — Spectral');

            propertyAnnotations(end+1) = CAFunctionProperty('rhoFunction', 'takes $$z$$ values and returns the no-motion density.');
            propertyAnnotations(end+1) = CAFunctionProperty('N2Function', 'takes $$z$$ values and returns the squared buoyancy frequency of the no-motion density.');
            propertyAnnotations(end+1) = CAPropertyAnnotation('verticalModes', 'instance of the InternalModes class');
            propertyAnnotations(end+1) = CANumericProperty('rho_nm0',{'z'},'kg m^{-3}', '$$\rho_\textrm{nm}(z)$$, no-motion density at time `t0`');
            propertyAnnotations(end+1) = CANumericProperty('N2',{'z'},'rad^2 s^{-2}', '$$N^2(z)$$, squared buoyancy frequency of the no-motion density, $$N^2\equiv - \frac{g}{\rho_0} \frac{\partial \rho_\textrm{nm}}{\partial z}$$');
            propertyAnnotations(end+1) = CANumericProperty('dLnN2',{'z'},'', '$$\frac{\partial \ln N^2}{\partial z}$$, vertical variation of the log of the squared buoyancy frequency');
            propertyAnnotations(end+1) = CANumericProperty('FinvMatrix',{'z','j'},'', 'transformation matrix $$F_g^{-1}$$');
            propertyAnnotations(end+1) = CANumericProperty('FMatrix',{'j','z'},'', 'transformation matrix $$F_g$$');
            propertyAnnotations(end+1) = CANumericProperty('GinvMatrix',{'z','j'},'', 'transformation matrix $$G_g^{-1}$$');
            propertyAnnotations(end+1) = CANumericProperty('GMatrix',{'j','z'},'', 'transformation matrix $$G_g$$');
            propertyAnnotations(end+1) = CANumericProperty('z_int',{'z'},'', 'Quadrature weights for the vertical grid');
            propertyAnnotations(end+1) = CANumericProperty('rho0',{},'kg m^{-3}', 'density of $$\rho_\textrm{nm}$$ at the surface (z=0)', detailedDescription='- topic: Domain Attributes');
            propertyAnnotations(end).attributes('standard_name') = 'sea_surface_density';

            propertyAnnotations(end+1) = CANumericProperty('Lz',{},'m', 'domain size in the z-direction', detailedDescription='- topic: Domain Attributes — Grid — Spatial');
            propertyAnnotations(end+1) = CANumericProperty('Nz',{},'', 'points in the z-coordinate, `length(z)`', detailedDescription='- topic: Domain Attributes — Grid — Spatial');
        end

        function [Lz,Nz,options] = requiredPropertiesForStratificationFromGroup(group,opts)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
                opts.shouldIgnoreMissingProperties logical = false
            end
            arguments (Output)
                Lz (1,1) double {mustBePositive}
                Nz (1,1) double {mustBePositive}
                options
            end
            vars = CAAnnotatedClass.propertyValuesFromGroup(group,WVStratification.namesOfRequiredPropertiesForStratification,shouldIgnoreMissingProperties=opts.shouldIgnoreMissingProperties);

            Nz = length(vars.z);
            Lz = vars.Lz;
            vars = rmfield(vars,{'Lz'});
            options = namedargs2cell(vars);
        end

        % function methodAnnotations = methodAnnotationsForStratification()
        %     % return array of WVAnnotations to annotate the methods
        %     %
        %     % This function returns annotations for all methods of the
        %     % WVStratification class.
        %     %
        %     % - Topic: Developer
        %     % - Declaration: methodAnnotations = WVStratification.methodAnnotationsForStratifiedFlow()
        %     % - Returns methodAnnotations: array of WVAnnotations instances
        %     methodAnnotations = WVAnnotation.empty(0,0);
        % 
        %     methodAnnotations(end+1) = WVAnnotation('diffZF', 'differentiates a variable of (x,y,z) by projecting onto the F-modes, differentiating, and transforming back to (x,y,z)');
        %     methodAnnotations(end+1) = WVAnnotation('diffZG', 'differentiates a variable of (x,y,z) by projecting onto the G-modes, differentiating, and transforming back to (x,y,z)');
        % end
    end

end