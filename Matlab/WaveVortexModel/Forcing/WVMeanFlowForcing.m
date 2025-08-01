classdef WVMeanFlowForcing < WVForcing
    % Resonant forcing at the natural frequency of each mode
    %
    % The unforced model basically looks likes like this,
    %
    % $$
    % \frac{\partial}{\partial t} A^{klj} = F_\textrm{inertial}^{klj} + F_\textrm{damp}^{klj}
    % $$
    %
    % for each of the three components. The forcing adds a new term,
    %
    % $$
    % \frac{\partial}{\partial t} A^{klj} = \underbrace{M_{A}^{klj} \left(\bar{A}^{klj}  - A^{klj} \right)/ \tau}_{F_\textrm{force}} + F_\textrm{inertial}^{klj} + F_\textrm{damp}^{klj}
    % $$
    %
    % which forces those select modes to relax to their $$\bar{A}^{klj}$$
    % state with time scale $$\tau$$.  If the time scale is set to 0, then the mean
    % amplitudes remain fixed for all time. In that limit, the
    % equations can be written as,
    %
    % $$
    % \frac{\partial}{\partial t} A^{klj} = \neg M_{A}^{klj} \left( F_\textrm{inertial}^{klj} + F_\textrm{damp}^{klj} \right)
    % $$
    %
    % This is most often used when initializing a model, e.g.,
    %
    % ```matlab
    % model = WVModel(wvt,nonlinearFlux=WVNonlinearFluxForced(wvt,uv_damp=wvt.uvMax));
    % ```
    %
    % - Topic: Initializing
    % - Declaration: WVNonlinearFluxForced < [WVNonlinearFlux](/classes/wvnonlinearflux/)
    properties
        A0_indices (:,1) uint64 = []    % Forcing mask, A0. 1s at the forced modes, 0s at the unforced modes
        Ap_indices (:,1) uint64 = []    % Forcing mask, Ap. 1s at the forced modes, 0s at the unforced modes
        Am_indices (:,1) uint64 = []    % Forcing mask, Am. 1s at the forced modes, 0s at the unforced modes

        A0bar (:,1) double = []  % A0 'mean' value to relax to
        Apbar (:,1) double = []  % Ap 'mean' value to relax to
        Ambar (:,1) double = []  % Am 'mean' value to relax to

        tau0 (1,1) double = 0    % A0 relaxation time
        tauP (1,1) double = 0    % Ap relaxation time
        tauM (1,1) double = 0    % Am relaxation time
    end

    methods
        function self = WVMeanFlowForcing(wvt,options)
            % initialize the WVNonlinearFlux nonlinear flux
            %
            % - Declaration: nlFlux = WVNonlinearFlux(wvt,options)
            % - Parameter wvt: a WVTransform instance
            % - Parameter uv_damp: (optional) characteristic speed used to set the damping. Try using wvt.uvMax.
            % - Parameter w_damp: (optional) characteristic speed used to set the damping. Try using wvt.wMax.
            % - Parameter nu_xy: (optional) coefficient for damping
            % - Parameter nu_z: (optional) coefficient for damping
            % - Returns nlFlux: a WVNonlinearFlux instance
            arguments
                wvt WVTransform {mustBeNonempty}
                options.name {mustBeText}

                options.Apbar (:,1) double = []
                options.Ambar (:,1) double = []
                options.A0bar (:,1) double = []
                options.A0_indices (:,1) uint64 = []
                options.Ap_indices (:,1) uint64 = []
                options.Am_indices (:,1) uint64 = []
                options.tauP (1,1) double = 0
                options.tauM (1,1) double = 0
                options.tau0 (1,1) double = 0
            end
            self@WVForcing(wvt,options.name,WVForcingType(["Spectral","PVSpectral"]));

            if ~isfield(options,"name")
                error("You must specify a unique name for these spectral masks, e.g., geostrophic mean flow, or M2 tide.")
            end

            canInitializeDirectly = any(options.A0_indices(:)) | any(options.Ap_indices(:)) | any(options.Am_indices(:));
            if canInitializeDirectly == true
                self.Apbar= options.Apbar;
                self.Ambar= options.Ambar;
                self.A0bar= options.A0bar;
                self.A0_indices  = options.A0_indices;
                self.Ap_indices  = options.Ap_indices;
                self.Am_indices  = options.Am_indices;
                self.tauP = options.tauP;
                self.tauM = options.tauM;
                self.tau0 = options.tau0;
            end
        end
        function setWaveForcingCoefficients(self,Apbar,Ambar,options)
            % set forcing values for the wave part of the flow
            %
            % Forcing takes the following form,
            %
            % $$
            % \frac{\partial}{\partial t} A_\pm^{klj} = \underbrace{M_{A_\pm}^{klj} \left(\bar{A}_\pm^{klj}  - A_\pm^{klj} \right)/ \tau_\pm}_{F_\textrm{force}} + F_\textrm{inertial}^{klj} + F_\textrm{damp}^{klj}
            % $$
            %
            % where $$M_{A_\pm}^{klj}$$ are masks (1s and 0s),
            % $$\bar{A}_\pm^{klj}$$ are mean amplitudes, and $$\tau_\pm$$
            % are time scales. If the time scale is set to 0, then the mean
            % amplitudes remain fixed for all time. In that limit, the
            % equations can be written as,
            %
            % $$
            % \frac{\partial}{\partial t} A_\pm^{klj} = \neg M_{A_\pm}^{klj} \left( F_\textrm{inertial}^{klj} + F_\textrm{damp}^{klj} \right)
            % $$
            %
            % - Topic: Set forcing
            % - Declaration:  setWaveForcingCoefficients(Apbar,Ambar,options)
            % - Parameter Apbar: Ap 'mean' value to relax to
            % - Parameter Ambar: Am 'mean' value to relax to
            % - Parameter MAp: (optional) forcing mask, Ap. 1s at the forced modes, 0s at the unforced modes. Default is MAp = abs(Apbar) > 0
            % - Parameter MAm: (optional) forcing mask, Am. 1s at the forced modes, 0s at the unforced modes. Default is MAm = abs(Apbar) > 0
            % - Parameter tauP: (optional) relaxation time (default 0)
            % - Parameter tauM: (optional) relaxation time (default 0)
            arguments
                self WVMeanFlowForcing {mustBeNonempty}
                Apbar (:,:) double {mustBeNonempty}
                Ambar (:,:) double {mustBeNonempty}
                options.MAp (:,:) logical = abs(Apbar) > 1e-6*max(abs(Apbar(:)))
                options.MAm (:,:) logical = abs(Ambar) > 1e-6*max(abs(Ambar(:)))
                options.tauP (1,1) double = 0
                options.tauM (1,1) double = 0
            end
            if self.wvt.hasForcingWithName("spectral vanishing viscosity")
                svv = self.wvt.forcingWithName("spectral vanishing viscosity");
                dampedIndicesAp = options.MAp(self.wvt.Kh > svv.k_damp);
                dampedIndicesAm = options.MAm(self.wvt.Kh > svv.k_damp);
                if any(dampedIndicesAp(:) | dampedIndicesAm(:))
                    warning('You have set %d forcing modes in the damping region. These will be removed.',sum(dampedIndicesAp(:))+sum(dampedIndicesAm(:)));
                    Apbar(self.wvt.Kh > svv.k_damp) = 0;
                    options.MAp(self.wvt.Kh > svv.k_damp) = 0;
                    Ambar(self.wvt.Kh > svv.k_damp) = 0;
                    options.MAm(self.wvt.Kh > svv.k_damp) = 0;
                end
            end

            self.Ap_indices = find(options.MAp);
            self.Apbar = Apbar(self.Ap_indices);
            self.tauP = options.tauP;

            self.Am_indices = find(options.MAm);
            self.Ambar = Ambar(self.Am_indices);
            self.tauM = options.tauM;

            fprintf('You are forcing at %d wave modes.\n',length(self.Ap_indices) + length(self.Am_indices));
        end

        function setGeostrophicForcingCoefficients(self,A0bar,options)
            % set forcing values for the geostrophic part of the flow
            %
            % Forcing takes the following form,
            %
            % $$
            % \frac{\partial}{\partial t} A_0^{klj} = \underbrace{M_{A_0}^{klj} \left(\bar{A}_0^{klj}  - A_0^{klj} \right)/ \tau_0}_{F_\textrm{force}} + F_\textrm{inertial}^{klj} + F_\textrm{damp}^{klj}
            % $$
            %
            % where $$M_{A_0}^{klj}$$ are masks (1s and 0s),
            % $$\bar{A}_0^{klj}$$ are mean amplitudes, and $$\tau_0$$
            % are time scales. If the time scale is set to 0, then the mean
            % amplitudes remain fixed for all time. In that limit, the
            % equations can be written as,
            %
            % $$
            % \frac{\partial}{\partial t} A_0^{klj} = \neg M_{A_0}^{klj} \left( F_\textrm{inertial}^{klj} + F_\textrm{damp}^{klj} \right)
            % $$
            %
            % - Topic: Set forcing
            % - Declaration: setGeostrophicForcingCoefficients(A0bar,options)
            % - Parameter A0bar: A0 'mean' value to relax to
            % - Parameter MA0: (optional) forcing mask, A0. 1s at the forced modes, 0s at the unforced modes. If it is left blank, then it will be produced using the nonzero values of A0bar
            % - Parameter tau0: (optional) relaxation time
            arguments
                self WVMeanFlowForcing {mustBeNonempty}
                A0bar (:,:) double {mustBeNonempty}
                options.MA0 (:,:) logical = abs(A0bar) > 1e-6*max(abs(A0bar(:)))
                options.tau0 (1,1) double = 0
            end

            if self.wvt.hasForcingWithName("spectral vanishing viscosity")
                svv = self.wvt.forcingWithName("spectral vanishing viscosity");
                dampedIndicesA0 = options.MA0(self.wvt.Kh > svv.k_damp);

                if any(dampedIndicesA0(:))
                    warning('You have set %d forcing modes in the damping region. These will be removed.',sum(dampedIndicesA0(:)));
                    A0bar(self.wvt.Kh > svv.k_damp) = 0;
                    options.MA0(self.wvt.Kh > svv.k_damp) = 0;
                end
            end

            self.A0_indices = find(options.MA0);
            self.A0bar = A0bar(self.A0_indices);
            self.tau0 = options.tau0;

            fprintf('You are forcing at %d geostrophic modes.\n',length(self.A0_indices));
        end

        function [model_spectrum, r] = setNarrowBandGeostrophicForcing(self, options)
            arguments
                self WVMeanFlowForcing {mustBeNonempty}
                options.r (1,1) double
                options.k_r (1,1) double =(self.wvt.k(2)-self.wvt.k(1))*2
                options.k_f (1,1) double =(self.wvt.k(2)-self.wvt.k(1))*20
                options.j_f (1,1) double = 1
                options.u_rms (1,1) double = 0.2 % set the *total* energy (not just kinetic) equal to 0.5*u_rms^2
                options.initialPV {mustBeMember(options.initialPV,{'none','narrow-band','full-spectrum'})} = 'narrow-band'
            end

            if ~isa(self.wvt,"WVGeometryDoublyPeriodicBarotropic")
                % the idea is to set the energy at the sea-surface and
                % so we need to know the relative amplitude of this
                % mode at the surface.
                F = self.wvt.FinvMatrix;
                surfaceMag = 1/F(end,options.j_f+1);
                sbRatio = abs(F(end,options.j_f+1)/F(1,options.j_f+1));
                % sbRatio = 1; % should we change the damping scale? Or no?
                h = self.wvt.h(options.j_f+1);
                magicNumber = 2.25;
            else
                surfaceMag = 1;
                sbRatio = 1;
                h = self.wvt.h;
                magicNumber = 0.0225;
            end


            if isfield(options,"r")
                r = options.r;
                k_r = self.r/(magicNumber*options.u_rms);
            else
                r = magicNumber*sbRatio*options.u_rms*options.k_r; % 1/s bracket [0.02 0.025]
                % fprintf('1/r is %.1f days, switching to %.1f days\n',1/(self.r*86400),1/(r*86400));
                k_r = options.k_r;
            end
            k_f = options.k_f;
            j_f = options.j_f;
            wvt = self.wvt;

            % smallDampIndex = find(abs(self.damp(:,1)) > 1.1*abs(self.r),1,'first');
            % fprintf('(k_r=%.2f dk, k_f=%d dk, k_nu=%d dk.\n',k_r/wvt.dk,round(k_f/wvt.dk),round(self.k_damp/wvt.dk));
            % fprintf('Small scale damping begins around k=%d dk. You have k_f=%d dk.\n',smallDampIndex-1,round(k_f/(wvt.k(2)-wvt.k(1))));


            deltaK = wvt.kRadial(2)-wvt.kRadial(1);
            MA0 = zeros(wvt.spectralMatrixSize);
            MA0(wvt.Kh > k_f-deltaK/2 & wvt.Kh < k_f+deltaK/2 & wvt.J == j_f) = 1;

            if strcmp(options.initialPV,'narrow-band') || strcmp(options.initialPV,'full-spectrum')
                u_rms = surfaceMag * options.u_rms;

                m = 3/2; % We don't really know what this number is.
                kappa_epsilon = 0.5 * u_rms^2 / ( ((3*m+5)/(2*m+2))*k_r^(-2/3) - k_f^(-2/3) );
                model_viscous = @(k) kappa_epsilon * k_r^(-5/3 - m) * k.^m;
                model_inverse = @(k) kappa_epsilon * k.^(-5/3);
                model_forward = @(k) kappa_epsilon * k_f^(4/3) * k.^(-3);
                model_spectrum = @(k) model_viscous(k) .* (k<k_r) + model_inverse(k) .* (k >= k_r & k<=k_f) + model_forward(k) .* (k>k_f);

                [~,~,wvt.A0] = wvt.geostrophicComponent.randomAmplitudesWithSpectrum(A0Spectrum= @(k,j) model_spectrum(k),shouldOnlyRandomizeOrientations=1);

                if strcmp(options.initialPV,'narrow-band')
                    wvt.A0 = MA0 .* wvt.A0;
                else
                    if isa(self.wvt,"WVGeometryDoublyPeriodicBarotropic")
                        u = wvt.u;
                        v = wvt.v;
                    else
                        u = wvt.ssu;
                        v = wvt.ssv;
                    end
                    zeta = wvt.ssh;
                    KE = mean(mean(0.5*(u.^2+v.^2)));
                    PE = mean(mean(0.5*(9.81*zeta.^2)/h));
                    u_rms_surface = mean(mean(sqrt(u.^2+v.^2)));
                    fprintf("surface u_rms: %.2g cm/s\n",100*u_rms_surface);
                    fprintf("surface energy, %g.\n",KE+PE);
                    fprintf('desired energy: %g, actual energy %g\n',0.5 * u_rms^2,wvt.geostrophicEnergy/h);
                end
            end
            self.setGeostrophicForcingCoefficients(MA0 .* wvt.A0,MA0=MA0,tau0=0);
        end
        
        function [Fp, Fm, F0] = addSpectralForcing(self, wvt, Fp, Fm, F0)
            if self.tauP > 0
                Fp(self.Ap_indices) = (self.Apbar - wvt.Ap(self.Ap_indices))/self.tauP + Fp(self.Ap_indices);
            elseif ~isempty(self.Ap_indices)
                Fp(self.Ap_indices) = 0;
                if ~isequal(wvt.Ap(self.Ap_indices),self.Apbar)
                    wvt.Ap(self.Ap_indices) = self.Apbar;
                end
            end
            if self.tauM > 0
                Fm(self.Am_indices) = (self.Ambar - wvt.Am(self.Am_indices))/self.tauM + Fm(self.Am_indices);
            elseif ~isempty(self.Am_indices)
                Fm(self.Am_indices) = 0;
                if ~isequal(wvt.Am(self.Am_indices),self.Ambar)
                    wvt.Am(self.Am_indices) = self.Ambar;
                end
            end
            if self.tau0 > 0
                F0(self.A0_indices) = (self.A0bar - wvt.A0(self.A0_indices))/self.tau0 + F0(self.A0_indices);
            elseif ~isempty(self.A0_indices)
                F0(self.A0_indices) = 0;
                if ~isequal(wvt.A0(self.A0_indices),self.A0bar)
                    wvt.A0(self.A0_indices) = self.A0bar;
                end
            end
        end

        function F0 = addPotentialVorticitySpectralForcing(self, wvt, F0)
            if self.tau0 > 0
                F0(self.A0_indices) = (self.A0bar - wvt.A0(self.A0_indices))/self.tau0 + F0(self.A0_indices);
            else
                F0(self.A0_indices) = 0;
                % This if-statement is a performance gain, given the
                % tradeoffs, and it will ensure the flux is doing as
                % expected.
                if ~isequal(wvt.A0(self.A0_indices),self.A0bar)
                    wvt.A0(self.A0_indices) = self.A0bar;
                end
            end
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            options.name = self.name;
            
            options.tauP = self.tauP;
            options.tauM = self.tauM;
            options.tau0 = self.tau0;

            if ~isempty(self.Ap_indices)
                Abar = zeros(self.wvt.spectralMatrixSize);
                Abar(self.Ap_indices) = self.Apbar;
                [AbarX2] = self.wvt.spectralVariableWithResolution(wvtX2,Abar);
                options.Apbar = self.Apbar;
                options.Ap_indices = find(AbarX2);
            end

            if ~isempty(self.Am_indices)
                Abar = zeros(self.wvt.spectralMatrixSize);
                Abar(self.Am_indices) = self.Ambar;
                [AbarX2] = self.wvt.spectralVariableWithResolution(wvtX2,Abar);
                options.Ambar = self.Ambar;
                options.Am_indices = find(AbarX2);
            end

            if ~isempty(self.A0_indices)
                Abar = zeros(self.wvt.spectralMatrixSize);
                Abar(self.A0_indices) = self.A0bar;
                [AbarX2] = self.wvt.spectralVariableWithResolution(wvtX2,Abar);
                options.A0bar = self.A0bar;
                options.A0_indices = find(AbarX2);
            end

            optionArgs = namedargs2cell(options);
            force = WVMeanFlowForcing(wvtX2,optionArgs{:});
        end

        function writeToGroup(self,group,propertyAnnotations,attributes)
            arguments
                self CAAnnotatedClass
                group NetCDFGroup
                propertyAnnotations CAPropertyAnnotation = CAPropertyAnnotation.empty(0,0)
                attributes = configureDictionary("string","string")
            end
            % override the logic, and only pass non-zero coefficients.
            properties = {'name'};
            if ~isempty(self.Ap_indices)
                properties = union(properties,{'Ap_indices','Apbar','tauP'});
            end
            if ~isempty(self.Am_indices)
                properties = union(properties,{'Am_indices','Ambar','tauM'});
            end
            if ~isempty(self.A0_indices)
                properties = union(properties,{'A0_indices','A0bar','tau0'});
            end
            propertyAnnotations = self.propertyAnnotationWithName(properties);      
            writeToGroup@CAAnnotatedClass(self,group,propertyAnnotations,attributes);
        end
    end

    methods (Static)
        function vars = classRequiredPropertyNames()
            vars = {'name','Ap_indices','Apbar','tauP','Am_indices','Ambar','tauM','A0_indices','A0bar','tau0'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CAPropertyAnnotation('name','name of the forcing');
            propertyAnnotations(end+1) = CADimensionProperty('Ap_indices', '','indices into the Ap matrix');
            propertyAnnotations(end+1) = CANumericProperty('Apbar', {'Ap_indices'}, '','Ap mean value',isComplex=true);
            propertyAnnotations(end+1) = CANumericProperty('tauP', {}, 's','Ap relaxation time');
            propertyAnnotations(end+1) = CADimensionProperty('Am_indices', '','indices into the Am matrix');
            propertyAnnotations(end+1) = CANumericProperty('Ambar', {'Am_indices'}, '','Am mean value',isComplex=true);
            propertyAnnotations(end+1) = CANumericProperty('tauM', {}, 's','Am relaxation time');
            propertyAnnotations(end+1) = CADimensionProperty('A0_indices', '','indices into the A0 matrix');
            propertyAnnotations(end+1) = CANumericProperty('A0bar', {'A0_indices'}, '','A0 mean value',isComplex=true);
            propertyAnnotations(end+1) = CANumericProperty('tau0', {}, 's','A0 relaxation time');
        end
    end

end