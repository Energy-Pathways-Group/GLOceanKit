classdef WVAdaptiveDamping < WVForcing
    % Adaptive small-scale damping
    %
    % This damping operator is a linear closure designed not to mix
    % geostrophic and wavemodes. This requires setting the diffusivity
    % equal to the viscosity. It also uses the spectral vanishing viscosity
    % filter to prevent any damping below a particular wavenumber.
    %
    %
    % - Topic: Initializing
    % - Declaration: WVAdaptiveDamping < [WVForcing](/classes/forcing/wvforcing/)
    properties
        damp
        k_damp % wavenumber at which the significant scale damping starts.
        k_no_damp % wavenumber below which there is zero damping
        j_damp % wavenumber at which the significant scale damping starts.
        j_no_damp % wavenumber below which there is zero damping

        forcingListener
        assumedEffectiveHorizontalGridResolution = Inf;
    end

    methods
        function self = WVAdaptiveDamping(wvt)
            % initialize the WVNonlinearFlux nonlinear flux
            %
            % Note that you should never need to set
            % shouldAssumeAntialiasing to true, because the WVGeometry will
            % hand back the correct effective grid resolution whether you
            % enabled it at the transform level, or as a filter.
            %
            % - Declaration: nlFlux = WVAdaptiveViscosity(wvt)
            % - Parameter wvt: a WVTransform instance
            % - Returns self: a WVAdaptiveViscosity instance
            arguments
                wvt WVTransform {mustBeNonempty}
            end
            self@WVForcing(wvt,"adaptive damping",WVForcingType(["Spectral","PVSpectral"]));
            self.wvt = wvt;
            self.isClosure = true;
            self.buildDampingOperator();
            self.forcingListener = addlistener(self.wvt,'forcingDidChange',@self.forcingDidChangeNotification);
        end

        function didGetRemovedFromTransform(self, wvt)
            delete(self.forcingListener);
            self.forcingListener = [];
        end

        function forcingDidChangeNotification(self,~,~)
            if self.wvt.effectiveHorizontalGridResolution ~= self.assumedEffectiveHorizontalGridResolution
                self.buildDampingOperator();
            end
        end

        function buildDampingOperator(self)
            % Builds the damping operator
            %
            % - Declaration: buildDampingOperator(self)
            % - Parameter self: an instance of WVAdaptiveViscosity
            % - Returns: None
            arguments
                self WVAdaptiveDamping {mustBeNonempty}
            end
            self.assumedEffectiveHorizontalGridResolution = self.wvt.effectiveHorizontalGridResolution;

            kl_max = pi/self.assumedEffectiveHorizontalGridResolution;
            j_max = self.wvt.effectiveJMax;
            j_index = find(self.wvt.j == self.wvt.effectiveJMax);
            [K,L,~] = self.wvt.kljGrid;
            [Qkl,Qj,self.k_no_damp,self.k_damp,self.j_no_damp,self.j_damp] = self.spectralVanishingViscosityFilter(kl_max, j_max);
            prefactor_xy = self.assumedEffectiveHorizontalGridResolution/(pi^2);
            prefactor_z = (pi*pi*self.wvt.Lr2(j_index)/(self.assumedEffectiveHorizontalGridResolution)^2)*prefactor_xy;

            Lr2inv = 1./self.wvt.Lr2;
            self.damp = -prefactor_xy*Qkl.*(K.^2 +L.^2) ;
            if ~isa(self.wvt,"WVGeometryDoublyPeriodicBarotropic")
                self.damp = self.damp - prefactor_z*Qj.*Lr2inv;
            end
        end

        function [Qkl,Qj,kl_cutoff, kl_damp, j_cutoff, j_damp] = spectralVanishingViscosityFilter(self, kl_max, j_max)
            % Builds the spectral vanishing viscosity operator
            %
            % - Declaration: spectralVanishingViscosityFilter(self, options)
            % - Parameter self: an instance of WVAdaptiveViscosity
            % - Parameter options: struct with field shouldAssumeAntialiasing
            % - Returns: Qkl, Qj, kl_cutoff, kl_damp
            arguments
                self WVAdaptiveDamping {mustBeNonempty}
                kl_max
                j_max
            end
            wvt_ = self.wvt;
            dkl_min = min(wvt_.dk, wvt_.dl);
            kl_cutoff = dkl_min*(kl_max/dkl_min)^(3/4);

            b = sqrt(-log(0.1));
            kl_damp = (kl_max+b*kl_cutoff)/(1+b); % approximately

            [K,L,J] = wvt_.kljGrid;
            Kh = sqrt(K.^2 + L.^2);

            Qkl = exp( - ((abs(Kh)-kl_max)./(abs(Kh)-kl_cutoff)).^2 );
            Qkl(abs(Kh)<kl_cutoff) = 0;
            Qkl(abs(Kh)>kl_max) = 1;

            if wvt_.Nj > 2
                dj = wvt_.j(2)-wvt_.j(1);
                j_cutoff = dj*(j_max/dj)^(3/4);
                j_damp = (j_max+b*j_cutoff)/(1+b); % approximately
                Qj = exp( - ((J-j_max)./(J-j_cutoff)).^2 );
                Qj(J<j_cutoff) = 0;
                Qj(J>j_max) = 1;
            else
                j_cutoff = 0;
                j_damp = 0;
                Qj = ones(size(J));
            end
        end

        % function [Qkl,Qj,kl_cutoff, kl_damp, j_cutoff, j_damp] = spectralVanishingViscosityFilter(self, options)
        %     % Builds the spectral vanishing viscosity operator
        %     %
        %     % - Declaration: spectralVanishingViscosityFilter(self, options)
        %     % - Parameter self: an instance of WVAdaptiveViscosity
        %     % - Parameter options: struct with field shouldAssumeAntialiasing
        %     % - Returns: Qkl, Qj, kl_cutoff, kl_damp
        %     arguments
        %         self WVAdaptiveDamping {mustBeNonempty}
        %         options.shouldAssumeAntialiasing logical = false
        %     end
        %     wvt_ = self.wvt;
        %     k_max = max(wvt_.k);
        %     l_max = max(wvt_.l);
        %     j_max = max(wvt_.j);
        %     if options.shouldAssumeAntialiasing == 1
        %         k_max = 2*k_max/3;
        %         l_max = 2*l_max/3;
        %         j_max = 2*j_max/3;
        %     end
        % 
        %     kl_max = min(k_max,l_max);
        %     dkl_min = min(wvt_.dk, wvt_.dl);
        %     kl_cutoff = dkl_min*(kl_max/dkl_min)^(3/4);
        % 
        %     b = sqrt(-log(0.1));
        %     kl_damp = (kl_max+b*kl_cutoff)/(1+b); % approximately
        % 
        %     [K,L,J] = wvt_.kljGrid;
        %     Kh = sqrt(K.^2 + L.^2);
        % 
        %     Qkl = exp( - ((abs(Kh)-kl_max)./(abs(Kh)-kl_cutoff)).^2 );
        %     Qkl(abs(Kh)<kl_cutoff) = 0;
        %     Qkl(abs(Kh)>kl_max) = 1;
        % 
        %     if wvt_.Nj > 2
        %         dj = wvt_.j(2)-wvt_.j(1);
        %         j_cutoff = dj*(j_max/dj)^(3/4);
        %         j_damp = (j_max+b*j_cutoff)/(1+b); % approximately
        %         Qj = exp( - ((J-j_max)./(J-j_cutoff)).^2 );
        %         Qj(J<j_cutoff) = 0;
        %         Qj(J>j_max) = 1;
        %     else
        %         j_cutoff = 0;
        %         j_damp = 0;
        %         Qj = ones(size(J));
        %     end
        % end

        function dampingTimeScale = dampingTimeScale(self)
            % Computes the damping time scale
            %
            % - Declaration: dampingTimeScale(self)
            % - Parameter self: an instance of WVAdaptiveViscosity
            % - Returns: dampingTimeScale
            arguments
                self WVAdaptiveDamping {mustBeNonempty}
            end
            dampingTimeScale = 1/max(abs(self.F0_damp(:)));
        end
        
        function [Fp, Fm, F0] = addSpectralForcing(self, wvt, Fp, Fm, F0)
            % Adds spectral forcing
            %
            % - Declaration: addSpectralForcing(self, wvt, Fp, Fm, F0)
            % - Parameter self: an instance of WVAdaptiveViscosity
            % - Parameter wvt: a WVTransform instance
            % - Parameter Fp: positive frequency forcing
            % - Parameter Fm: negative frequency forcing
            % - Parameter F0: zero frequency forcing
            % - Returns: Fp, Fm, F0
            arguments
                self WVAdaptiveDamping {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
                Fp double {mustBeNonempty}
                Fm double {mustBeNonempty}
                F0 double {mustBeNonempty}
            end
            uvMax = wvt.uvMax;
            Fp = Fp + uvMax * self.damp .* wvt.Ap;
            Fm = Fm + uvMax * self.damp .* wvt.Am;
            F0 = F0 + uvMax * self.damp .* wvt.A0;
        end

        function F0 = addPotentialVorticitySpectralForcing(self, wvt, F0)
            % Adds potential vorticity spectral forcing
            %
            % - Declaration: addPotentialVorticitySpectralForcing(self, wvt, F0)
            % - Parameter self: an instance of WVAdaptiveViscosity
            % - Parameter wvt: a WVTransform instance
            % - Parameter F0: zero frequency forcing
            % - Returns: F0
            arguments
                self WVAdaptiveDamping {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
                F0 double {mustBeNonempty}
            end
            F0 = F0 + wvt.uvMax * self.damp .* wvt.A0;
        end

        function force = forcingWithResolutionOfTransform(self, wvtX2)
            % Creates a forcing with the resolution of the transform
            %
            % - Declaration: forcingWithResolutionOfTransform(self, wvtX2)
            % - Parameter self: an instance of WVAdaptiveViscosity
            % - Parameter wvtX2: a WVTransform instance with doubled resolution
            % - Returns: force
            arguments
                self WVAdaptiveDamping {mustBeNonempty}
                wvtX2 WVTransform {mustBeNonempty}
            end
            force = WVAdaptiveDamping(wvtX2);
        end
    end
    methods (Static)
        function vars = classRequiredPropertyNames()
            % Returns the required property names for the class
            %
            % - Declaration: classRequiredPropertyNames()
            % - Returns: vars
            arguments
            end
            vars = {};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            % Returns the defined property annotations for the class
            %
            % - Declaration: classDefinedPropertyAnnotations()
            % - Returns: propertyAnnotations
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
        end
    end
end