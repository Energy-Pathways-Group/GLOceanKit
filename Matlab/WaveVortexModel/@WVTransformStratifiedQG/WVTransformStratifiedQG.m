classdef WVTransformStratifiedQG < WVGeometryDoublyPeriodicStratified & WVTransform & WVGeostrophicMethods
    % A class for disentangling hydrostatic waves and vortices in variable stratification
    %
    % To initialization an instance of the WVTransformHydrostatic class you
    % must specific the domain size, the number of grid points and *either*
    % the density profile or the stratification profile.
    %
    % ```matlab
    % N0 = 3*2*pi/3600;
    % L_gm = 1300;
    % N2 = @(z) N0*N0*exp(2*z/L_gm);
    % wvt = WVTransformHydrostatic([100e3, 100e3, 4000],[64, 64, 65], N2=N2,latitude=30);
    % ```
    %
    % - Topic: Initialization
    % - Topic: Primary flow components
    % - Topic: Stratification
    % - Topic: Stratification — Vertical modes
    % - Topic: Stratification — Validation
    % - Topic: Initial conditions
    % - Topic: Energetics of flow components
    % - Topic: Operations
    %
    % - Declaration: classdef WVTransformHydrostatic < [WVTransform](/classes/wvtransform/)
    properties (Dependent)
        totalEnergySpatiallyIntegrated
        totalEnergy
        isHydrostatic
    end

    properties (GetAccess=private,SetAccess=private)
        Fpv, F0
    end

    methods
        function self = WVTransformStratifiedQG(Lxyz, Nxyz, options)
            % create a wave-vortex transform for variable stratification
            %
            % Creates a new instance of the WVTransformHydrostatic class
            % appropriate for disentangling hydrostatic waves and vortices
            % in variable stratification
            %
            % You must initialization by passing *either* the density
            % profile or the stratification profile.
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
                Lxyz (1,3) double {mustBePositive}
                Nxyz (1,3) double {mustBePositive}
                options.shouldAntialias (1,1) logical = true
                options.z (:,1) double {mustBeNonempty} % quadrature points!
                options.j (:,1) double {mustBeNonempty}
                options.Nj (1,1) double {mustBePositive}
                options.rhoFunction function_handle = @isempty
                options.N2Function function_handle = @isempty
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.planetaryRadius (1,1) double = 6.371e6
                options.rotationRate (1,1) double = 7.2921E-5
                options.latitude (1,1) double = 33
                options.g (1,1) double = 9.81
                
                options.dLnN2 (:,1) double
                options.PF0inv
                options.QG0inv
                options.PF0
                options.QG0
                options.h_0 (:,1) double
                options.P0 (:,1) double
                options.Q0 (:,1) double
                options.z_int (:,1) double
            end

            optionArgs = namedargs2cell(options);
            self@WVGeometryDoublyPeriodicStratified(Lxyz, Nxyz, optionArgs{:})
            self@WVTransform(WVForcingType(["PVSpectral","PVSpatial","PVSpectralAmplitude"]));
            self@WVGeostrophicMethods();

            self.initializeGeostrophicComponent();

            % This is not good, I think this should go in the constructor.
            self.addForcing(WVNonlinearAdvection(self));

            % the property annotations for these variables will already
            % have beena added, but that is okay, they will be replaced.
            varNames = self.namesOfTransformVariables();
            self.addOperation(self.operationForKnownVariable(varNames{:}),shouldOverwriteExisting=true,shouldSuppressWarning=true);

            self.A0 = zeros(self.spectralMatrixSize);
            self.F0 = zeros(self.spectralMatrixSize);
            self.Fpv = zeros(self.spatialMatrixSize);
            if self.geostrophicComponent.normalization ~= "qgpv"
                error("This transform requires the geostrophic component to be normalized the the qgpv norm.");
                % self.A0PV = self.geostrophicComponent.multiplierForVariable(WVCoefficientMatrix.A0,"qgpv-inv");
            end
        end

        function wvtX2 = waveVortexTransformWithResolution(self,m)
            names = {'shouldAntialias','N2Function','rho0','planetaryRadius','rotationRate','latitude','g'};
            optionArgs = {};
            for i=1:length(names)
                optionArgs{2*i-1} = names{i};
                optionArgs{2*i} = self.(names{i});
            end
            wvtX2 = WVTransformStratifiedQG([self.Lx self.Ly self.Lz],m,optionArgs{:});
            forcing = WVForcing.empty(0,length(self.forcing));
            for iForce=1:length(self.forcing)
                forcing(iForce) = self.forcing(iForce).forcingWithResolutionOfTransform(wvtX2);
            end
            wvtX2.setForcing(forcing);

            wvtX2.t0 = self.t0;
            wvtX2.t = self.t;
            [wvtX2.A0] = self.spectralVariableWithResolution(wvtX2,self.A0);
        end

        function wvt = hydrostaticTransform(self)
            names = {'shouldAntialias','N2Function','rho0','planetaryRadius','rotationRate','latitude','g'};
            optionArgs = {};
            for i=1:length(names)
                optionArgs{2*i-1} = names{i};
                optionArgs{2*i} = self.(names{i});
            end
            wvt = WVTransformHydrostatic([self.Lx self.Ly self.Lz],[self.Nx self.Ny self.Nz],optionArgs{:});
            forcing = WVForcing.empty(0,length(self.forcing));
            for iForce=1:length(self.forcing)
                forcing(iForce) = self.forcing(iForce).forcingWithResolutionOfTransform(wvt);
            end
            wvt.setForcing(forcing);

            wvt.t0 = self.t0;
            wvt.t = self.t;
            wvt.A0 = self.spectralVariableWithResolution(wvt,self.A0);
        end

        function energy = get.totalEnergySpatiallyIntegrated(self)
            [u,v,eta] = self.variableWithName('u','v','eta');
            energy = sum(shiftdim(self.z_int,-2).*mean(mean( u.^2 + v.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
        end

        function energy = get.totalEnergy(self)
            energy = sum( self.A0_TE_factor(:).*( abs(self.A0(:)).^2) );
        end
        
        function flag = get.isHydrostatic(self)
            flag = true;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Nonlinear flux computation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function A0 = transformQGPVToWaveVortex(self,qgpv)
            A0 = self.transformFromSpatialDomainWithFg(self.transformFromSpatialDomainWithFourier(qgpv));
        end

        function F0 = nonlinearFlux(self)
            self.Fpv = 0*self.Fpv;
            for i=1:length(self.spatialFluxForcing)
                self.Fpv = self.spatialFluxForcing(i).addPotentialVorticitySpatialForcing(self,self.Fpv);
            end
            self.F0 = self.transformQGPVToWaveVortex(self.Fpv);
            for i=1:length(self.spectralFluxForcing)
                self.F0 = self.spectralFluxForcing(i).addPotentialVorticitySpectralForcing(self,self.F0);
            end
            for i=1:length(self.spectralAmplitudeForcing)
                self.F0 = self.spectralAmplitudeForcing(i).setPotentialVorticitySpectralForcing(self,self.F0);
            end
            F0 = self.F0;
        end

        function F0 = fluxForForcing(self)
            arguments (Input)
                self WVTransform
            end
            arguments (Output)
                F0 dictionary
            end
            F0 = configureDictionary("string","cell");
            self.Fpv = 0*self.Fpv;
            for i=1:length(self.spatialFluxForcing)
               Fpv0 = self.Fpv;
               self.Fpv = self.spatialFluxForcing(i).addPotentialVorticitySpatialForcing(self,self.Fpv);
               F0{self.spatialFluxForcing(i).name} = self.transformQGPVToWaveVortex(self.Fpv-Fpv0);
            end
            self.F0 = self.transformQGPVToWaveVortex(self.Fpv);
            for i=1:length(self.spectralFluxForcing)
                F0_i = self.F0;
                self.F0 = self.spectralFluxForcing(i).addPotentialVorticitySpectralForcing(self,self.F0);
                F0{self.spectralFluxForcing(i).name} = self.F0 - F0_i;
            end
            for i=1:length(self.spectralAmplitudeForcing)
                F0_i = self.F0;
                self.F0 = self.spectralAmplitudeForcing(i).setPotentialVorticitySpectralForcing(self,self.F0);
                F0{self.spectralAmplitudeForcing(i).name} = self.F0 - F0_i;
            end
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
            u = self.transformToSpatialDomainWithFourier(self.PF0inv*(self.P0 .* options.A0));
        end

        function w = transformToSpatialDomainWithG(self, options)
            arguments
                self WVTransform {mustBeNonempty}
                options.Apm double = 0
                options.A0 double = 0
            end
            w = self.transformToSpatialDomainWithFourier(self.QG0inv*(self.Q0 .* options.A0));
        end

    end

    methods (Static)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % CAAnnotatedClass required methods, which enables writeToFile
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = WVTransformStratifiedQG.propertyAnnotationsForTransform();
        end

        function vars = classRequiredPropertyNames()
            vars = WVTransformStratifiedQG.namesOfRequiredPropertiesForTransform();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Stratification specific property annotations and initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function requiredPropertyNames = namesOfRequiredPropertiesForTransform()
            requiredPropertyNames = WVGeometryDoublyPeriodicStratified.namesOfRequiredPropertiesForGeometry();
            requiredPropertyNames = union(requiredPropertyNames,WVTransformStratifiedQG.newRequiredPropertyNames);
        end

        function newRequiredPropertyNames = newRequiredPropertyNames()
            newRequiredPropertyNames = {'A0','kl','t0','t','forcing'};
        end

        function names = namesOfTransformVariables()
            names = {'A0t','uvMax','zeta_z','ssh','ssu','ssv','u','v','eta','pi','p','psi','qgpv','rho_e','rho_total'};
        end

        function propertyAnnotations = propertyAnnotationsForTransform()
            spectralDimensionNames = WVTransformStratifiedQG.spectralDimensionNames();
            spatialDimensionNames = WVTransformStratifiedQG.spatialDimensionNames();

            propertyAnnotations = WVGeometryDoublyPeriodicStratified.propertyAnnotationsForGeometry();
            propertyAnnotations = cat(2,propertyAnnotations,WVGeostrophicMethods.propertyAnnotationsForGeostrophicComponent(spectralDimensionNames = spectralDimensionNames));
            transformProperties = WVTransform.propertyAnnotationsForTransform('A0','A0_TE_factor','A0_QGPV_factor','A0_TZ_factor',spectralDimensionNames = spectralDimensionNames);

            varNames = WVTransformStratifiedQG.namesOfTransformVariables();
            varAnnotations = WVTransform.propertyAnnotationForKnownVariable(varNames{:},spectralDimensionNames = spectralDimensionNames,spatialDimensionNames = spatialDimensionNames);
            propertyAnnotations = cat(2,propertyAnnotations,transformProperties,varAnnotations);
        end

        function [Lxyz,Nxyz,options] = requiredPropertiesForTransformFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                Lxyz (1,3) double {mustBePositive}
                Nxyz (1,3) double {mustBePositive}
                options
            end
            [Lxyz, Nxyz, geomOptions] = WVGeometryDoublyPeriodicStratified.requiredPropertiesForGeometryFromGroup(group);
            % CAAnnotatedClass.throwErrorIfMissingProperties(group,WVTransformBarotropicQG.newRequiredPropertyNames);
            % vars = CAAnnotatedClass.propertyValuesFromGroup(group,WVTransformBarotropicQG.newRequiredPropertyNames);
            % newOptions = namedargs2cell(vars);
            % options = cat(2,geomOptions,newOptions);
            options = geomOptions;
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
                options.shouldReadOnly logical = false
            end
            arguments (Output)
                wvt WVTransform
                ncfile NetCDFFile
            end
            ncfile = NetCDFFile(path,shouldReadOnly=options.shouldReadOnly);
            wvt = WVTransformStratifiedQG.transformFromGroup(ncfile);
            wvt.initFromNetCDFFile(ncfile,iTime=options.iTime,shouldDisplayInit=1);
            wvt.initForcingFromNetCDFFile(ncfile);
        end


        function wvt = transformFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                wvt WVTransform {mustBeNonempty}
            end  
            [Lxy, Nxy, options] = WVTransformStratifiedQG.requiredPropertiesForTransformFromGroup(group);
            wvt = WVTransformStratifiedQG(Lxy,Nxy,options{:});
        end

    end

end



