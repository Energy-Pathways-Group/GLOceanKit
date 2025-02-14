classdef WVSpectralMasks < WVForcing
    % 3D forced nonlinear flux for Boussinesq flow
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
        MA0 = []    % Forcing mask, A0. 1s at the forced modes, 0s at the unforced modes
        MAp = []    % Forcing mask, Ap. 1s at the forced modes, 0s at the unforced modes
        MAm = []    % Forcing mask, Am. 1s at the forced modes, 0s at the unforced modes

        A0bar = []  % A0 'mean' value to relax to
        Apbar = []  % Ap 'mean' value to relax to
        Ambar = []  % Am 'mean' value to relax to

        tau0 = 0    % A0 relaxation time
        tauP = 0    % Ap relaxation time
        tauM = 0    % Am relaxation time
    end

    methods
        function self = WVSpectralMasks(wvt,options)
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
                options.name {mustBeText} = "spectral masks"

                options.Apbar (:,:) double
                options.Ambar (:,:) double
                options.A0bar (:,:) double
                options.MAp   (:,:) logical
                options.MAm   (:,:) logical
                options.MA0   (:,:) logical 
                options.tauP (1,1) double
                options.tauM (1,1) double
                options.tau0 (1,1) double
            end
            self@WVForcing(wvt,options.name,WVForcingType(["Spectral","PVSpectral"]));

            canInitializeDirectly = all(isfield(options, self.classRequiredPropertyNames));
            if canInitializeDirectly == true
                self.Apbar= options.Apbar;
                self.Ambar= options.Ambar;
                self.A0bar= options.A0bar;
                self.MAp  = options.MAp;
                self.MAm  = options.MAm;
                self.MA0  = options.MA0;
                self.tauP = options.tauP;
                self.tauM = options.tauM;
                self.tau0 = options.tau0;
            else
                self.Apbar= zeros(wvt.spectralMatrixSize);
                self.Ambar= zeros(wvt.spectralMatrixSize);
                self.A0bar= zeros(wvt.spectralMatrixSize);
                self.MAp  = false(wvt.spectralMatrixSize);
                self.MAm  = false(wvt.spectralMatrixSize);
                self.MA0  = false(wvt.spectralMatrixSize);
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
                self WVSpectralMasks {mustBeNonempty}
                Apbar (:,:) double {mustBeNonempty}
                Ambar (:,:) double {mustBeNonempty}
                options.MAp (:,:) logical = abs(Apbar) > 0
                options.MAm (:,:) logical = abs(Ambar) > 0
                options.tauP (1,1) double = 0
                options.tauM (1,1) double = 0
            end
            dampedIndicesAp = options.MAp(self.wvt.Kh > self.k_damp);
            dampedIndicesAm = options.MAm(self.wvt.Kh > self.k_damp);
            warning('You have set %d forcing modes in the damping region. These will be removed.',sum(dampedIndicesAp(:))+sum(dampedIndicesAm(:)));
            Apbar(self.wvt.Kh > self.k_damp) = 0;
            options.MAp(self.wvt.Kh > self.k_damp) = 0;
            Ambar(self.wvt.Kh > self.k_damp) = 0;
            options.MAm(self.wvt.Kh > self.k_damp) = 0;

            self.Apbar = Apbar;
            self.MAp = options.MAp;
            self.tauP = options.tauP;

            self.Ambar = Ambar;
            self.MAm = options.MAm;
            self.tauM = options.tauM;

            fprintf('You are forcing at %d wave modes.\n',sum(options.MAp(:))+sum(options.MAm(:)));
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
                self WVSpectralMasks {mustBeNonempty}
                A0bar (:,:) double {mustBeNonempty}
                options.MA0 (:,:) logical = abs(A0bar) > 0
                options.tau0 (1,1) double = 0
            end

            self.A0bar = A0bar;
            self.MA0 = options.MA0;
            self.tau0 = options.tau0;
        end
        
        function [Fp, Fm, F0] = addSpectralForcing(self, wvt, Fp, Fm, F0)
            if self.tauP > 0
                Fp = self.MAp.*(self.Apbar - wvt.Ap)/self.tauP + Fp;
            elseif ~isempty(self.MAp)
                Fp = (~self.MAp) .* Fp;
            end
            if self.tauM > 0
                Fm = self.MAm.*(self.Ambar - wvt.Am)/self.tauM + Fm;
            elseif ~isempty(self.MAm)
                Fm = (~self.MAm) .* Fm;
            end
            if self.tau0 > 0
                F0 = self.MA0.*(self.A0bar - wvt.A0)/self.tau0 + F0;
            elseif ~isempty(self.MA0)
                F0 = (~self.MA0) .* F0;
            end
        end

        function F0 = addPotentialVorticitySpectralForcing(self, wvt, F0)
            if self.tau0 > 0
                F0 = self.MA0.*(self.A0bar - wvt.A0)/self.tau0 + F0;
            elseif ~isempty(self.MA0)
                F0 = (~self.MA0) .* F0;
            end
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            options.tauP = self.tauP;
            options.tauM = self.tauM;
            options.tau0 = self.tau0;
            [options.MAp, options.Apbar] = self.wvt.spectralVariableWithResolution(wvtX2,self.MAp, self.Apbar);
            [options.MAm, options.Ambar] = self.wvt.spectralVariableWithResolution(wvtX2,self.MAm, self.Ambar);
            [options.MA0, options.A0bar] = self.wvt.spectralVariableWithResolution(wvtX2,self.MA0, self.A0bar);
            optionArgs = namedargs2cell(options);
            force = WVSpectralMasks(wvtX2,optionArgs{:});
        end
    end

    methods (Static)
        function vars = classRequiredPropertyNames()
            vars = {'name','MAp','Apbar','tauP','MAm','Ambar','tauM','MA0','A0bar','tau0'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CAPropertyAnnotation('name','name of the forcing');
            propertyAnnotations(end+1) = CANumericProperty('MAp', {'j','kl'}, '','Ap mask');
            propertyAnnotations(end+1) = CANumericProperty('Apbar', {'j','kl'}, '','Ap mean value');
            propertyAnnotations(end+1) = CANumericProperty('tauP', {}, 's','Ap relaxation time');
            propertyAnnotations(end+1) = CANumericProperty('MAm', {'j','kl'}, '','Am mask');
            propertyAnnotations(end+1) = CANumericProperty('Ambar', {'j','kl'}, '','Am mean value');
            propertyAnnotations(end+1) = CANumericProperty('tauM', {}, 's','Am relaxation time');
            propertyAnnotations(end+1) = CANumericProperty('MA0', {'j','kl'}, '','A0 mask');
            propertyAnnotations(end+1) = CANumericProperty('A0bar', {'j','kl'}, '','A0 mean value');
            propertyAnnotations(end+1) = CANumericProperty('tau0', {}, 's','A0 relaxation time');
        end
    end

end