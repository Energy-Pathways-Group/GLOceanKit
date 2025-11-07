classdef WVTransformConstantStratification < WVGeometryDoublyPeriodicStratifiedConstant & WVTransform & WVGeostrophicMethods & WVInternalGravityWaveMethods & WVInertialOscillationMethods & WVMeanDensityAnomalyMethods
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
    end
    properties
        Fu, Fv, Feta
        nonlinearFluxFunction

        cos_alpha
        sin_alpha
        ApmD_scaled
        ApmW_scaled
    end

    methods
        function self = WVTransformConstantStratification(Lxyz, Nxyz, options)
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

                options.N0 (1,1) double {mustBePositive} = 5.2e-3
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.planetaryRadius (1,1) double = 6.371e6
                options.rotationRate (1,1) double = 7.2921E-5
                options.latitude (1,1) double = 33
                options.g (1,1) double = 9.81
                
                options.isHydrostatic logical = false
            end

            optionArgs = namedargs2cell(options);
            self@WVGeometryDoublyPeriodicStratifiedConstant(Lxyz, Nxyz, optionArgs{:})
            self@WVTransform(WVForcingType(["HydrostaticSpatial","Spectral","SpectralAmplitude"]));
            self@WVGeostrophicMethods();
            self@WVMeanDensityAnomalyMethods();
            self@WVInternalGravityWaveMethods();
            self@WVInertialOscillationMethods();

            self.initializeGeostrophicComponent();
            self.initializeMeanDensityAnomalyComponent();
            self.initializeInternalGravityWaveComponent();
            self.initializeInertialOscillationComponent();

            % This is not good, I think this should go in the constructor.
            self.addForcing(WVNonlinearAdvection(self));

            % the property annotations for these variables will already
            % have beena added, but that is okay, they will be replaced.
            varNames = self.namesOfTransformVariables();
            self.addOperation(self.operationForKnownVariable(varNames{:}),shouldOverwriteExisting=true,shouldSuppressWarning=true);

            self.addOperation(self.operationForKnownVariable('u','v','w','eta','p',flowComponent=self.geostrophicComponent));
            self.addOperation(self.operationForKnownVariable('u','v','w','eta','p',flowComponent=self.waveComponent));
            self.addOperation(self.operationForKnownVariable('u','v','w','eta','p',flowComponent=self.inertialComponent));
            self.addOperation(self.operationForKnownVariable('u','v','w','eta','p',flowComponent=self.mdaComponent));
            self.addOperation(EtaTrueOperation());
            self.addOperation(APVOperation());
            self.addOperation(APEOperation(self));
            
            self.A0 = zeros(self.spectralMatrixSize);
            self.Ap = zeros(self.spectralMatrixSize);
            self.Am = zeros(self.spectralMatrixSize);

            self.Fu=zeros(self.spatialMatrixSize);
            self.Fv=zeros(self.spatialMatrixSize);
            self.Feta=zeros(self.spatialMatrixSize);

            k = shiftdim(self.k,-1);
            l = shiftdim(self.l,-1);
            kappa = sqrt(k.^2 + l.^2);
            self.cos_alpha = k./kappa;
            self.sin_alpha = l./kappa;
            self.cos_alpha(1) = 0;
            self.sin_alpha(1) = 0;

            signNorm = -2*(mod(self.j,2) == 1)+1; % equivalent to (-1)^j
            prefactor = signNorm * sqrt((self.g*self.Lz)/(2*(self.N0*self.N0 - self.f*self.f)));
            mj = (self.j*pi/self.Lz);
            self.ApmD_scaled = (mj/2) .* prefactor;
            self.ApmW_scaled = sqrt(-1) * (kappa/2) .* prefactor;

            if self.isHydrostatic
                self.nonlinearFluxFunction = @() self.nonlinearFluxHydrostatic();
            else
                self.nonlinearFluxFunction = @() self.nonlinearFluxNonhydrostatic();
            end
        end

        function wvtX2 = waveVortexTransformWithResolution(self,m)
            names = {'shouldAntialias','N2Function','rho0','planetaryRadius','rotationRate','latitude','g'};
            optionArgs = {};
            for i=1:length(names)
                optionArgs{2*i-1} = names{i};
                optionArgs{2*i} = self.(names{i});
            end
            wvtX2 = WVTransformConstantStratification([self.Lx self.Ly self.Lz],m,optionArgs{:});
            forcing = WVForcing.empty(0,length(self.forcing));
            for iForce=1:length(self.forcing)
                forcing(iForce) = self.forcing(iForce).forcingWithResolutionOfTransform(wvtX2);
            end
            wvtX2.setForcing(forcing);

            wvtX2.t0 = self.t0;
            wvtX2.t = self.t;
            [wvtX2.A0,wvtX2.Ap,wvtX2.Am] = self.spectralVariableWithResolution(wvtX2,self.A0,self.Ap,self.Am);
        end

        function wvt2 = waveVortexTransformWithExplicitAntialiasing(self)
            if self.shouldAntialias == false
                error("This function only applies to transforms that are dealiasing.")
            end
            names = {'shouldAntialias','N2Function','rho0','planetaryRadius','rotationRate','latitude','g'};
            optionArgs = {};
            for i=1:length(names)
                optionArgs{2*i-1} = names{i};
                optionArgs{2*i} = self.(names{i});
                if names{i} == "shouldAntialias"
                    optionArgs{2*i} = false;
                end
            end
            wvt2 = WVTransformConstantStratification([self.Lx self.Ly self.Lz],[self.Nx self.Ny self.Nz],optionArgs{:});
            wvt2.removeAllForcing();
            wvt2.addForcing(WVAntialiasing(wvt2));

            for iForce=1:length(self.forcing)
                wvt2.addForcing(self.forcing(iForce).forcingWithResolutionOfTransform(wvt2));
            end

            wvt2.t0 = self.t0;
            wvt2.t = self.t;
            [wvt2.A0,wvt2.Ap,wvt2.Am] = self.spectralVariableWithResolution(wvt2,self.A0,self.Ap,self.Am);

        end

        function dx = effectiveHorizontalGridResolution(self)
            %returns the effective grid resolution in meters
            %
            % The effective grid resolution is the highest fully resolved
            % wavelength in the model. This value takes into account
            % anti-aliasing, and is thus appropriate for setting damping
            % operators.
            %
            % - Topic: Properties
            % - Declaration: flag = effectiveHorizontalGridResolution(other)
            % - Returns effectiveHorizontalGridResolution: double
            arguments
                self WVGeometryDoublyPeriodic
            end
            if self.hasForcingWithName("antialias filter")
                dx = self.forcingWithName("antialias filter").effectiveHorizontalGridResolution;
            else
                dx = effectiveHorizontalGridResolution@WVGeometryDoublyPeriodic(self);
            end
        end

        function j_max = effectiveJMax(self)
            if self.hasForcingWithName("antialias filter")
                j_max = self.forcingWithName("antialias filter").effectiveJMax;
            else
                j_max =effectiveJMax@WVGeometryDoublyPeriodicStratified(self);
            end
        end

        function energy = get.totalEnergySpatiallyIntegrated(self)
            if self.isHydrostatic == 1
                [u,v,eta] = self.variableWithName('u','v','eta');
                energy = sum(shiftdim(self.z_int,-2).*mean(mean( u.^2 + v.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
            else
                [u,v,w,eta] = self.variableWithName('u','v','w','eta');
                energy = sum(shiftdim(self.z_int,-2).*mean(mean( u.^2 + v.^2 + w.^2 + shiftdim(self.N2,-2).*eta.*eta, 1 ),2 ) )/2;
            end
        end

        function energy = get.totalEnergy(self)
            energy = sum( self.Apm_TE_factor(:).*( abs(self.Ap(:)).^2 + abs(self.Am(:)).^2 ) + self.A0_TE_factor(:).*( abs(self.A0(:)).^2) );
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Nonlinear flux computation
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % function [Fp,Fm,F0] = nonlinearFlux(self)
        %     % self.Fu(:)=0;self.Fv(:)=0;self.Feta(:)=0;
        %     % self.Fu=0*self.Fu;self.Fv(:)=0*self.Fv;self.Feta(:)=0*self.Feta;
        %     self.Fu=zeros(self.spatialMatrixSize);self.Fv=zeros(self.spatialMatrixSize);self.Feta=zeros(self.spatialMatrixSize);
        %     for i=1:length(self.spatialFluxForcing)
        %        [self.Fu, self.Fv, self.Feta] = self.spatialFluxForcing(i).addHydrostaticSpatialForcing(self, self.Fu, self.Fv, self.Feta);
        %     end
        %     [Fp,Fm,F0] = self.transformUVEtaToWaveVortex(self.Fu, self.Fv, self.Feta);
        %     for i=1:length(self.spectralFluxForcing)
        %        [Fp,Fm,F0] = self.spectralFluxForcing(i).addSpectralForcing(self,Fp, Fm, F0);
        %     end
        %     for i=1:length(self.spectralAmplitudeForcing)
        %        [Fp,Fm,F0] = self.spectralAmplitudeForcing(i).setSpectralForcing(self,Fp, Fm, F0);
        %     end  
        % end
        function [Fp,Fm,F0] = nonlinearFlux(self)
            [Fp,Fm,F0] = self.nonlinearFluxFunction();
        end
        function [Fp,Fm,F0] = nonlinearFluxHydrostatic(self)
            Fu=zeros(self.spatialMatrixSize);Fv=zeros(self.spatialMatrixSize);Feta=zeros(self.spatialMatrixSize); % this isn't good, need to cached
            for i=1:length(self.spatialFluxForcing)
                [Fu, Fv, Feta] = self.spatialFluxForcing(i).addHydrostaticSpatialForcing(self, Fu, Fv, Feta);
            end
            [Fp,Fm,F0] = self.transformUVEtaToWaveVortex(Fu, Fv, Feta);
            for i=1:length(self.spectralFluxForcing)
                [Fp,Fm,F0] = self.spectralFluxForcing(i).addSpectralForcing(self,Fp, Fm, F0);
            end
            for i=1:length(self.spectralAmplitudeForcing)
                [Fp,Fm,F0] = self.spectralAmplitudeForcing(i).setSpectralForcing(self,Fp, Fm, F0);
            end
        end

        function [Fp,Fm,F0] = nonlinearFluxNonhydrostatic(self)
            Fu=zeros(self.spatialMatrixSize);Fv=zeros(self.spatialMatrixSize);Fw=zeros(self.spatialMatrixSize);Feta=zeros(self.spatialMatrixSize); % this isn't good, need to cached
            for i=1:length(self.spatialFluxForcing)
                [Fu, Fv, Fw, Feta] = self.spatialFluxForcing(i).addNonhydrostaticSpatialForcing(self, Fu, Fv, Fw, Feta);
            end
            [Fp,Fm,F0] = self.transformUVWEtaToWaveVortex(Fu, Fv, Fw, Feta);
            for i=1:length(self.spectralFluxForcing)
                [Fp,Fm,F0] = self.spectralFluxForcing(i).addSpectralForcing(self,Fp, Fm, F0);
            end
            for i=1:length(self.spectralAmplitudeForcing)
                [Fp,Fm,F0] = self.spectralAmplitudeForcing(i).setSpectralForcing(self,Fp, Fm, F0);
            end
        end

        % function F = fluxForForcing(self)
        %     arguments (Input)
        %         self WVTransform
        %     end
        %     arguments (Output)
        %         F dictionary
        %     end
        %     F = configureDictionary("string","cell");
        %     Fu=0;Fv=0;Feta=0; % this isn't good, need to cached
        %     for i=1:length(self.spatialFluxForcing)
        %         Fu0=Fu;Fv0=Fv;Feta0=Feta;
        %         [Fu, Fv, Feta] = self.spatialFluxForcing(i).addHydrostaticSpatialForcing(self, Fu, Fv, Feta);
        %         [Fp,Fm,F0] = self.transformUVEtaToWaveVortex(Fu-Fu0, Fv-Fv0, Feta-Feta0);
        %         F{self.spatialFluxForcing(i).name} = struct("Fp",Fp,"Fm",Fm,"F0",F0);
        %     end
        %     [Fp,Fm,F0] = self.transformUVEtaToWaveVortex(Fu, Fv, Feta);
        %     for i=1:length(self.spectralFluxForcing)
        %         Fp_i = Fp; Fm_i = Fm; F0_i = F0;
        %         [Fp,Fm,F0] = self.spectralFluxForcing(i).addSpectralForcing(self,Fp, Fm, F0);
        %         F{self.spectralFluxForcing(i).name} = struct("Fp",Fp-Fp_i,"Fm",Fm-Fm_i,"F0",F0-F0_i);
        %     end
        %     for i=1:length(self.spectralAmplitudeForcing)
        %         Fp_i = Fp; Fm_i = Fm; F0_i = F0;
        %         [Fp,Fm,F0] = self.spectralAmplitudeForcing(i).setSpectralForcing(self,Fp, Fm, F0);
        %         F{self.spectralAmplitudeForcing(i).name} = struct("Fp",Fp-Fp_i,"Fm",Fm-Fm_i,"F0",F0-F0_i);
        %     end
        % end

        function varargout = spatialFluxForForcingWithName(self,name)
            Fu_name = replace(replace(join( ["Fu_", string(name)],"")," ","_"),"-","_");
            Fv_name = replace(replace(join( ["Fv_", string(name)],"")," ","_"),"-","_");
            Feta_name = replace(replace(join( ["Feta_", string(name)],"")," ","_"),"-","_");
            [Fu_,Fv_,Feta_] = self.variableWithName(Fu_name,Fv_name,Feta_name);
            if self.isHydrostatic && nargout == 3
                varargout = {Fu_,Fv_,Feta_};
            elseif self.isHydrostatic == false && nargout == 4
                Fw_name = replace(replace(join( ["Fw_", string(name)],"")," ","_"),"-","_");
                F_w = self.variableWithName(Fw_name);
                varargout = {Fu_,Fv_,F_w,Feta_};
            end
        end


    end

    methods (Static)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % CAAnnotatedClass required methods, which enables writeToFile
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = WVTransformConstantStratification.propertyAnnotationsForTransform();
        end

        function vars = classRequiredPropertyNames()
            vars = WVTransformConstantStratification.namesOfRequiredPropertiesForTransform();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Stratification specific property annotations and initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function requiredPropertyNames = namesOfRequiredPropertiesForTransform()
            requiredPropertyNames = WVGeometryDoublyPeriodicStratifiedConstant.namesOfRequiredPropertiesForGeometry();
            requiredPropertyNames = union(requiredPropertyNames,WVTransformConstantStratification.newRequiredPropertyNames);
        end

        function newRequiredPropertyNames = newRequiredPropertyNames()
            newRequiredPropertyNames = {'A0','Ap','Am','kl','t0','t','forcing'};
        end

        function names = namesOfTransformVariables()
            names = {'phase','conjPhase','A0t','Apt','Amt','uvMax','wMax','zeta_z','ssh','ssu','ssv','u','v','w','eta','pi','p','psi','qgpv','rho_e','rho_total'};
        end

        function propertyAnnotations = propertyAnnotationsForTransform()
            spectralDimensionNames = WVTransformConstantStratification.spectralDimensionNames();
            spatialDimensionNames = WVTransformConstantStratification.spatialDimensionNames();

            propertyAnnotations = WVGeometryDoublyPeriodicStratifiedConstant.propertyAnnotationsForGeometry();
            propertyAnnotations = cat(2,propertyAnnotations,WVGeostrophicMethods.propertyAnnotationsForGeostrophicComponent(spectralDimensionNames = spectralDimensionNames));
            transformProperties = WVTransform.propertyAnnotationsForTransform('A0','Ap','Am','A0_TE_factor','A0_QGPV_factor','A0_TZ_factor','A0_Psi_factor','Apm_TE_factor',spectralDimensionNames = spectralDimensionNames);

            varNames = WVTransformConstantStratification.namesOfTransformVariables();
            varAnnotations = WVTransform.propertyAnnotationForKnownVariable(varNames{:},spectralDimensionNames = spectralDimensionNames,spatialDimensionNames = spatialDimensionNames);
            propertyAnnotations = cat(2,propertyAnnotations,transformProperties,varAnnotations);
        end

      function [Lxyz, Nxyz, options] = requiredPropertiesForTransformFromGroup(group)
            arguments (Input)
                group NetCDFGroup {mustBeNonempty}
            end
            arguments (Output)
                Lxyz (1,3) double {mustBePositive}
                Nxyz (1,3) double {mustBePositive}
                options
            end
            [Lxyz, Nxyz, geomOptions] = WVGeometryDoublyPeriodicStratifiedConstant.requiredPropertiesForGeometryFromGroup(group);
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
            wvt = WVTransformConstantStratification.transformFromGroup(ncfile);
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
            [Lxy, Nxy, options] = WVTransformConstantStratification.requiredPropertiesForTransformFromGroup(group);
            wvt = WVTransformConstantStratification(Lxy,Nxy,options{:});
        end

    end

end



