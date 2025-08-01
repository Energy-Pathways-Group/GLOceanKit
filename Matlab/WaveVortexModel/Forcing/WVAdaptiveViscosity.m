classdef WVAdaptiveViscosity < WVForcing
    % Adaptive small-scale damping
    %
    % The damping is a simple Laplacian, but with a spectral vanishing
    % viscosity (SVV) operator applied that prevents any damping below a
    % cutoff wavenumber, adapted to the current fluid velocity. The SVV operator adjusts the wavenumbers being
    % damped depending on whether anti-aliasing is applied.
    %
    %
    % - Topic: Initializing
    % - Declaration: WVAdaptiveViscosity < [WVForcing](/classes/forcing/wvforcing/)
    properties
        Fpm_damp
        F0_damp

        k_damp % wavenumber at which the significant scale damping starts.
        k_no_damp % wavenumber below which there is zero damping
    end

    methods
        function self = WVAdaptiveViscosity(wvt,options)
            % initialize the WVNonlinearFlux nonlinear flux
            %
            % - Declaration: nlFlux = WVAdaptiveViscosity(wvt)
            % - Parameter wvt: a WVTransform instance
            % - Returns self: a WVAdaptiveViscosity instance
            arguments
                wvt WVTransform {mustBeNonempty}
                options.shouldAssumeAntialiasing logical = false
            end
            self@WVForcing(wvt,"adaptive svv",WVForcingType(["Spectral","PVSpectral"]));
            self.wvt = wvt;
            self.isClosure = true;
            self.buildDampingOperator(shouldAssumeAntialiasing=options.shouldAssumeAntialiasing);
        end

        function buildDampingOperator(self,options)
            % Builds the damping operator
            %
            % - Declaration: buildDampingOperator(self)
            % - Parameter self: an instance of WVAdaptiveViscosity
            % - Returns: None
            arguments
                self WVAdaptiveViscosity {mustBeNonempty}
                options.shouldAssumeAntialiasing logical = false
            end
            [K,L,J] = self.wvt.kljGrid;
            M = J*pi/self.wvt.Lz;
            [Qkl,Qj,self.k_no_damp,self.k_damp] = self.spectralVanishingViscosityFilter(shouldAssumeAntialiasing=options.shouldAssumeAntialiasing);
            prefactor_xy = self.wvt.effectiveHorizontalGridResolution/(pi^2);
            prefactor_z = self.wvt.effectiveVerticalGridResolution/(pi^2);
            if options.shouldAssumeAntialiasing == true
                prefactor_xy = 3*prefactor_xy/2;
                prefactor_z = 3*prefactor_z/2;
            end
            % (max(self.wvt.N2)/self.wvt.f/self.wvt.f)*
            Lr2inv = 1./self.wvt.Lr2;
            % self.damp = -prefactor_xy*(Qkl.*(K.^2 +L.^2) + Qj.*Lr2inv);
            if ~isempty(intersect(self.wvt.forcingType,WVForcingType("Spectral")))
                self.Fpm_damp = -prefactor_xy*Qkl.*(K.^2 +L.^2);
            end
            RV = self.wvt.geostrophicComponent.multiplierForVariable(WVCoefficientMatrix.A0,"rv");
            self.F0_damp = -prefactor_xy*Qkl.*RV.*(K.^2 +L.^2);
        end

        function [Qkl,Qj,kl_cutoff, kl_damp] = spectralVanishingViscosityFilter(self, options)
            % Builds the spectral vanishing viscosity operator
            %
            % - Declaration: spectralVanishingViscosityFilter(self, options)
            % - Parameter self: an instance of WVAdaptiveViscosity
            % - Parameter options: struct with field shouldAssumeAntialiasing
            % - Returns: Qkl, Qj, kl_cutoff, kl_damp
            arguments
                self WVAdaptiveViscosity {mustBeNonempty}
                options.shouldAssumeAntialiasing logical = false
            end
            wvt_ = self.wvt;
            k_max = max(wvt_.k);
            l_max = max(wvt_.l);
            j_max = max(wvt_.j);
            if options.shouldAssumeAntialiasing == 1
                k_max = 2*k_max/3;
                l_max = 2*l_max/3;
                j_max = 2*j_max/3;
            end

            kl_max = min(k_max,l_max);
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
                Qj = exp( - ((J-j_max)./(J-j_cutoff)).^2 );
                Qj(J<j_cutoff) = 0;
                Qj(J>j_max) = 1;
            else
                Qj = ones(size(J));
            end
        end

        function dampingTimeScale = dampingTimeScale(self)
            % Computes the damping time scale
            %
            % - Declaration: dampingTimeScale(self)
            % - Parameter self: an instance of WVAdaptiveViscosity
            % - Returns: dampingTimeScale
            arguments
                self WVAdaptiveViscosity {mustBeNonempty}
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
                self WVAdaptiveViscosity {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
                Fp double {mustBeNonempty}
                Fm double {mustBeNonempty}
                F0 double {mustBeNonempty}
            end
            uvMax = wvt.uvMax;
            Fp = Fp + uvMax * self.Fpm_damp .* wvt.Ap;
            Fm = Fm + uvMax * self.Fpm_damp .* wvt.Am;
            F0 = F0 + uvMax * self.F0_damp .* wvt.A0;
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
                self WVAdaptiveViscosity {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
                F0 double {mustBeNonempty}
            end
            F0 = F0 + wvt.uvMax * self.F0_damp .* wvt.A0;
        end

        function force = forcingWithResolutionOfTransform(self, wvtX2)
            % Creates a forcing with the resolution of the transform
            %
            % - Declaration: forcingWithResolutionOfTransform(self, wvtX2)
            % - Parameter self: an instance of WVAdaptiveViscosity
            % - Parameter wvtX2: a WVTransform instance with doubled resolution
            % - Returns: force
            arguments
                self WVAdaptiveViscosity {mustBeNonempty}
                wvtX2 WVTransform {mustBeNonempty}
            end
            force = WVAdaptiveViscosity(wvtX2);
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