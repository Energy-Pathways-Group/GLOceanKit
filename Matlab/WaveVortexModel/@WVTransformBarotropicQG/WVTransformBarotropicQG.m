classdef WVTransformBarotropicQG < WVGeometryDoublyPeriodicBarotropic & WVGeostrophicMethods
    % A transform for modeling single-layer quasigeostrophic flow
    %
    % This is a two-dimensional, single-layer which may be interpreted as
    % the sea-surface height. The 'h' parameter is the equivalent depth,
    % and 0.80 m is a typical value for the first baroclinic mode.
    %
    % ```matlab
    % Lxy = 50e3;
    % Nxy = 256;
    % latitude = 25;
    % wvt = WVTransformSingleMode([Lxy, Lxy], [Nxy, Nxy], h=0.8, latitude=latitude);
    % ```
    %
    % - Topic: Initialization
    %
    % - Declaration: classdef WVTransformBarotropicQG < [WVTransform](/classes/wvtransform/)
    properties (Dependent)
        h_0
        totalEnergySpatiallyIntegrated
        totalEnergy
    end

    properties (GetAccess=private,SetAccess=private)
        Fpv, F0, A0PV
    end

    methods
        function self = WVTransformBarotropicQG(Lxy, Nxy, options)
            % create geometry for 2D barotropic flow
            %
            % ```matlab
            % Lxy = 50e3;
            % Nxy = 256;
            % wvt = Cartesian2DBarotropic([Lxy, Lxy], [Nxy, Nxy]);
            % ```
            %
            % - Topic: Initialization
            % - Declaration: wvt = Cartesian2DBarotropic(Lxyz, Nxyz, options)
            % - Parameter Lxy: length of the domain (in meters) in the two coordinate directions, e.g. [Lx Ly]
            % - Parameter Nxy: number of grid points in the two coordinate directions, e.g. [Nx Ny]
            % - Parameter shouldAntialias: (optional) whether or not to de-alias for quadratic multiplications
            % - Returns wvt: a new Cartesian2DBarotropic instance
            arguments
                Lxy (1,2) double {mustBePositive}
                Nxy (1,2) double {mustBePositive}
                options.shouldAntialias (1,1) logical = true
                options.rotationRate (1,1) double = 7.2921E-5
                options.latitude (1,1) double = 33
                options.g (1,1) double = 9.81
                options.h (1,1) double = 0.8
                options.j (1,1) double {mustBeMember(options.j,[0 1])} = 1
            end
            optionCell = namedargs2cell(options);
            self@WVGeometryDoublyPeriodicBarotropic(Lxy,Nxy,optionCell{:});
            self@WVGeostrophicMethods();
            
            self.initializeGeostrophicComponent();

            % the property annotations for these variables will already
            % have beena added, but that is okay, they will be replaced.
            varNames = WVTransformBarotropicQG.namesOfTransformVariables();
            self.addOperation(self.operationForKnownVariable(varNames{:}),shouldOverwriteExisting=true);

            self.A0 = zeros(self.spectralMatrixSize);
            self.F0 = zeros(self.spectralMatrixSize);
            self.Fpv = zeros(self.spatialMatrixSize);
            self.A0PV = 1./self.A0_QGPV_factor;
            self.A0PV(isinf(self.A0PV))=0;
            %self.addOperation(self.operationForKnownVariable('u','v','eta','p',flowComponent=self.geostrophicComponent));
        end

        function val = get.h_0(self)
            val = self.h;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Nonlinear flux computation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function F0 = nonlinearFlux(self)
            self.Fpv = 0*self.Fpv;
            for i=1:length(self.spatialForcing)
               self.Fpv = self.spatialForcing(i).addPotentialVorticitySpatialForcing(self,self.Fpv);
            end
            self.F0 = self.A0PV .* self.transformFromSpatialDomainWithFourier(self.Fpv);
            for i=1:length(self.spectralForcing)
               self.F0 = self.spectralForcing(i).addPotentialVorticitySpectralForcing(self,self.F0);
            end
            F0 = self.F0;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations FROM the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function u = transformFromSpatialDomainWithFg(~, u)
        end

        function w = transformFromSpatialDomainWithGg(~, w)
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations TO the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function u = transformToSpatialDomainWithF(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            u = self.transformToSpatialDomainWithFourier(options.A0);
        end

        function w = transformToSpatialDomainWithG(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            w = self.transformToSpatialDomainWithFourier(options.A0);
        end

        function wvtX2 = waveVortexTransformWithResolution(self,m)
            wvtX2 = WVTransformBarotropicQG([self.Lx self.Ly],m,h=self.h,latitude=self.latitude,rotationRate=self.rotationRate,g=self.g);
            wvtX2.t0 = self.t0;
            wvtX2.t = self.t;
            [wvtX2.Ap,wvtX2.Am,wvtX2.A0] = self.spectralVariableWithResolution(wvtX2,self.Ap,self.Am,self.A0);
            % wvtX2.nonlinearFluxOperation = self.nonlinearFluxOperation.nonlinearFluxWithResolutionOfTransform(wvtX2);
        end

        function ratio = maxFg(self,k0, l0, j0)
            ratio = 1;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Energetics (total)
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function energy = get.totalEnergySpatiallyIntegrated(self)
            [u,v,eta] = self.variableWithName('u','v','eta');
            energy = sum(shiftdim(self.z_int,-2).*mean(mean( u.^2 + v.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
        end

        function energy = get.totalEnergy(self)
            energy = sum( self.A0_TE_factor(:).*( abs(self.A0(:)).^2) );
        end

        summarizeDegreesOfFreedom(self)

        function setSSH(self,ssh,options)
            arguments
                self WVTransformBarotropicQG
                ssh
                options.shouldRemoveMeanPressure double {mustBeMember(options.shouldRemoveMeanPressure,[0 1])} = 0
            end
            if options.shouldRemoveMeanPressure == 1
                sshbar = mean(mean(ssh(self.X,self.Y)));
            else
                sshbar = 0;
            end
            psi = @(X,Y,Z) (self.g/self.f)*(ssh(X,Y)-sshbar);

            self.setGeostrophicStreamfunction(psi);
        end
    end

    methods (Static)

        function names = classDefinedSpectralDimensionNames()
            % return a cell array of property names required by the class
            %
            % This function returns an array of property names required to be written
            % by the class, in order to restore its state.
            %
            % - Topic: Developer
            % - Declaration:  names = classRequiredPropertyNames()
            % - Returns names: array strings
            arguments (Output)
                names cell
            end
            names = {'kl'};
        end

        function names = classDefinedSpatialDimensionNames()
            % return a cell array of the spatial dimension names
            %
            % This function returns an array of dimension names
            %
            % - Topic: Developer
            % - Declaration:  names = classDefinedSpatialDimensionNames()
            % - Returns names: array strings
            arguments (Output)
                names cell
            end
            names = {'x','y'};
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % CAAnnotatedClass required methods, which enables writeToFile
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = WVTransformBarotropicQG.propertyAnnotationsForTransform();
        end

        function vars = classRequiredPropertyNames()
            vars = WVTransformBarotropicQG.namesOfRequiredPropertiesForTransform();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Stratification specific property annotations and initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function requiredPropertyNames = namesOfRequiredPropertiesForTransform()
            requiredPropertyNames = WVGeometryDoublyPeriodicBarotropic.namesOfRequiredPropertiesForGeometry();
            requiredPropertyNames = union(requiredPropertyNames,WVTransformBarotropicQG.newRequiredPropertyNames);
        end

        function newRequiredPropertyNames = newRequiredPropertyNames()
            newRequiredPropertyNames = {'A0','kl','t0','t'};
        end

        function names = namesOfTransformVariables()
            names = {'A0t','uvMax','zeta_z','ssh','u','v','eta','pi','psi','qgpv'};
        end

        function propertyAnnotations = propertyAnnotationsForTransform()
            spectralDimensionNames = WVTransformBarotropicQG.classDefinedSpectralDimensionNames();
            spatialDimensionNames = WVTransformBarotropicQG.classDefinedSpatialDimensionNames();

            propertyAnnotations = WVGeometryDoublyPeriodicBarotropic.propertyAnnotationsForGeometry();
            propertyAnnotations = cat(2,propertyAnnotations,WVGeostrophicMethods.propertyAnnotationsForGeostrophicComponent(spectralDimensionNames = spectralDimensionNames));
            transformProperties = WVTransform.propertyAnnotationsForTransform('A0','A0_TE_factor','A0_QGPV_factor','A0_TZ_factor',spectralDimensionNames = spectralDimensionNames);

            varNames = WVTransformBarotropicQG.namesOfTransformVariables();
            varAnnotations = WVTransform.propertyAnnotationForKnownVariable(varNames{:},spectralDimensionNames = spectralDimensionNames,spatialDimensionNames = spatialDimensionNames);
            propertyAnnotations = cat(2,propertyAnnotations,transformProperties,varAnnotations);
        end

        function [Lxy,Nxy,options] = requiredPropertiesForTransformFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                Lxy (1,2) double {mustBePositive}
                Nxy (1,2) double {mustBePositive}
                options
            end
            [Lxy, Nxy, geomOptions] = WVGeometryDoublyPeriodicBarotropic.requiredPropertiesForGeometryFromGroup(group);
            vars = CAAnnotatedClass.propertyValuesFromGroup(group,WVTransformBarotropicQG.newRequiredPropertyNames);
            newOptions = namedargs2cell(vars);
            options = cat(2,geomOptions,newOptions);
        end

        function [wvt,ncfile] = waveVortexTransformFromFile(path,options)
            % Initialize a WVTransformHydrostatic instance from an existing file
            %
            % This static method is called by WVTransform.waveVortexTransformFromFile
            % and should not need to be called directly.
            %
            % - Topic: Initialization (Static)
            % - Declaration: wvt = waveVortexTransformFromFile(path,options)
            % - Parameter path: path to a NetCDF file
            % - Parameter iTime: (optional) time index to initialize from (default 1)
            arguments (Input)
                path char {mustBeFile}
                options.iTime (1,1) double {mustBePositive} = 1
            end
            arguments (Output)
                wvt WVTransform
                ncfile NetCDFFile
            end
            ncfile = NetCDFFile(path);
            wvt = WVTransformBarotropicQG.transformFromGroup(ncfile);
            wvt.initFromNetCDFFile(ncfile,iTime=options.iTime,shouldDisplayInit=1);
        end


        function wvt = transformFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                wvt WVTransform {mustBeNonempty}
            end
            CAAnnotatedClass.throwErrorIfMissingProperties(group,WVTransformBarotropicQG.namesOfRequiredPropertiesForTransform);
            [Lxy, Nxy, options] = WVGeometryDoublyPeriodicBarotropic.requiredPropertiesForGeometryFromGroup(group);
            wvt = WVTransformBarotropicQG(Lxy,Nxy,options{:});
        end

    end
end