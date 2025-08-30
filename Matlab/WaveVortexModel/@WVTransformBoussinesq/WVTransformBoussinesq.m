classdef WVTransformBoussinesq < WVGeometryDoublyPeriodicStratifiedBoussinesq & WVTransform & WVGeostrophicMethods & WVInternalGravityWaveMethods & WVInertialOscillationMethods & WVMeanDensityAnomalyMethods
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
    properties
        Fu, Fv, Feta
        delta_uhat
        delta_vhat
        ApmW
        Ddelta
    end

    methods
        function self = WVTransformBoussinesq(Lxyz, Nxyz, options)
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
                options.PFpmInv
                options.QGpmInv
                options.PFpm
                options.QGpm
                options.h_pm
                options.Ppm
                options.Qpm
                options.QGwg
                options.K2unique
                options.iK2unique
            end

            optionArgs = namedargs2cell(options);
            self@WVGeometryDoublyPeriodicStratifiedBoussinesq(Lxyz, Nxyz, optionArgs{:})
            self@WVTransform(WVForcingType(["NonhydrostaticSpatial","Spectral","SpectralAmplitude"]));
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

            self.A0 = zeros(self.spectralMatrixSize);
            self.Ap = zeros(self.spectralMatrixSize);
            self.Am = zeros(self.spectralMatrixSize);

            self.Fu=zeros(self.spatialMatrixSize);
            self.Fv=zeros(self.spatialMatrixSize);
            self.Feta=zeros(self.spatialMatrixSize);

            % the following operators have dimensions (z,kl)
            k = shiftdim(self.k,-1);
            l = shiftdim(self.l,-1);
            kappa = sqrt(k.^2 + l.^2);
            self.delta_uhat = k./kappa/2;
            self.delta_vhat = l./kappa/2;
            self.delta_uhat(1) = 0;
            self.delta_vhat(1) = 0;
            self.ApmW = self.g*kappa./(self.N2 - self.f*self.f)/2;

            self.Ddelta = (self.N2./(self.N2 - self.f*self.f)) .* self.QG0inv*(squeeze(self.Q0 ./ self.P0).*self.PF0);  
        end

        function wvtX2 = waveVortexTransformWithResolution(self,m)
            names = {'shouldAntialias','N2Function','rho0','planetaryRadius','rotationRate','latitude','g'};
            optionArgs = {};
            for i=1:length(names)
                optionArgs{2*i-1} = names{i};
                optionArgs{2*i} = self.(names{i});
            end
            wvtX2 = WVTransformBoussinesq([self.Lx self.Ly self.Lz],m,optionArgs{:});
            forcing = WVForcing.empty(0,length(self.forcing));
            for iForce=1:length(self.forcing)
                forcing(iForce) = self.forcing(iForce).forcingWithResolutionOfTransform(wvtX2);
            end
            wvtX2.setForcing(forcing);

            wvtX2.t0 = self.t0;
            wvtX2.t = self.t;
            [wvtX2.A0,wvtX2.Ap,wvtX2.Am] = self.spectralVariableWithResolution(wvtX2,self.A0,self.Ap,self.Am);
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

        function flag = get.isHydrostatic(self)
            flag = false;
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

        function F = fluxForForcing(self)
            arguments (Input)
                self WVTransform
            end
            arguments (Output)
                F dictionary
            end
            F = configureDictionary("string","cell");
            Fu=0;Fv=0;Fw=0;Feta=0; % this isn't good, need to cached
            for i=1:length(self.spatialFluxForcing)
                Fu0=Fu;Fv0=Fv;Fw0=Fw;Feta0=Feta;
                [Fu, Fv, Fw, Feta] = self.spatialFluxForcing(i).addNonhydrostaticSpatialForcing(self, Fu, Fv, Fw, Feta);
                [Fp,Fm,F0] = self.transformUVWEtaToWaveVortex(Fu-Fu0, Fv-Fv0, Fw-Fw0, Feta-Feta0);
                F{self.spatialFluxForcing(i).name} = struct("Fp",Fp,"Fm",Fm,"F0",F0);
            end
            [Fp,Fm,F0] = self.transformUVWEtaToWaveVortex(Fu, Fv, Fw, Feta);
            for i=1:length(self.spectralFluxForcing)
                Fp_i = Fp; Fm_i = Fm; F0_i = F0;
                [Fp,Fm,F0] = self.spectralFluxForcing(i).addSpectralForcing(self,Fp, Fm, F0);
                F{self.spectralFluxForcing(i).name} = struct("Fp",Fp-Fp_i,"Fm",Fm-Fm_i,"F0",F0-F0_i);
            end
            for i=1:length(self.spectralAmplitudeForcing)
                Fp_i = Fp; Fm_i = Fm; F0_i = F0;
                [Fp,Fm,F0] = self.spectralAmplitudeForcing(i).setSpectralForcing(self,Fp, Fm, F0);
                F{self.spectralAmplitudeForcing(i).name} = struct("Fp",Fp-Fp_i,"Fm",Fm-Fm_i,"F0",F0-F0_i);
            end
        end

        function [Fu,Fv,Fw,Feta] = spatialFluxForForcingWithName(self,name)
            Fu_name = replace(replace(join( ["Fu_", string(name)],"")," ","_"),"-","_");
            Fv_name = replace(replace(join( ["Fv_", string(name)],"")," ","_"),"-","_");
            Fw_name = replace(replace(join( ["Fw_", string(name)],"")," ","_"),"-","_");
            Feta_name = replace(replace(join( ["Feta_", string(name)],"")," ","_"),"-","_");
            [Fu,Fv,Fw,Feta] = self.variableWithName(Fu_name,Fv_name,Fw_name,Feta_name);
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
            wvt2 = WVTransformBoussinesq([self.Lx self.Ly self.Lz],[self.Nx self.Ny self.Nz],optionArgs{:});
            wvt2.removeAllForcing();
            wvt2.addForcing(WVAntialiasing(wvt2));

            for iForce=1:length(self.forcing)
                wvt2.addForcing(self.forcing(iForce).forcingWithResolutionOfTransform(wvt2));
            end

            wvt2.t0 = self.t0;
            wvt2.t = self.t;
            [wvt2.A0,wvt2.Ap,wvt2.Am] = self.spectralVariableWithResolution(wvt2,self.A0,self.Ap,self.Am);

        end


    end

    methods (Static)

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % CAAnnotatedClass required methods, which enables writeToFile
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function propertyAnnotations = classDefinedPropertyAnnotations()
            propertyAnnotations = WVTransformBoussinesq.propertyAnnotationsForTransform();
        end

        function vars = classRequiredPropertyNames()
            vars = WVTransformBoussinesq.namesOfRequiredPropertiesForTransform();
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Stratification specific property annotations and initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function requiredPropertyNames = namesOfRequiredPropertiesForTransform()
            requiredPropertyNames = WVGeometryDoublyPeriodicStratifiedBoussinesq.namesOfRequiredPropertiesForGeometry();
            requiredPropertyNames = union(requiredPropertyNames,WVTransformBoussinesq.newRequiredPropertyNames);
        end

        function newRequiredPropertyNames = newRequiredPropertyNames()
            newRequiredPropertyNames = {'A0','Ap','Am','kl','t0','t','forcing'};
        end

        function names = namesOfTransformVariables()
            names = {'phase','conjPhase','A0t','Apt','Amt','uvMax','wMax','zeta_x','zeta_y','zeta_z','ssh','ssu','ssv','u','v','w','eta','pi','p','psi','qgpv','rho_e','rho_total'};
        end

        function propertyAnnotations = propertyAnnotationsForTransform()
            spectralDimensionNames = WVTransformBoussinesq.spectralDimensionNames();
            spatialDimensionNames = WVTransformBoussinesq.spatialDimensionNames();

            propertyAnnotations = WVGeometryDoublyPeriodicStratifiedBoussinesq.propertyAnnotationsForGeometry();
            propertyAnnotations = cat(2,propertyAnnotations,WVGeostrophicMethods.propertyAnnotationsForGeostrophicComponent(spectralDimensionNames = spectralDimensionNames));
            transformProperties = WVTransform.propertyAnnotationsForTransform('A0','Ap','Am','A0_TE_factor','A0_QGPV_factor','A0_TZ_factor','A0_Psi_factor','Apm_TE_factor',spectralDimensionNames = spectralDimensionNames);

            varNames = WVTransformBoussinesq.namesOfTransformVariables();
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
            [Lxyz, Nxyz, geomOptions] = WVGeometryDoublyPeriodicStratifiedBoussinesq.requiredPropertiesForGeometryFromGroup(group);
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
            end
            arguments (Output)
                wvt WVTransform
                ncfile NetCDFFile
            end
            ncfile = NetCDFFile(path);
            wvt = WVTransformBoussinesq.transformFromGroup(ncfile);
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
            [Lxy, Nxy, options] = WVTransformBoussinesq.requiredPropertiesForTransformFromGroup(group);
            wvt = WVTransformBoussinesq(Lxy,Nxy,options{:});
        end

    end

end



