classdef WVAdaptiveViscosity < WVForcing
    % Small-scale damping
    %
    % The damping is a simple Laplacian, but with a spectral vanishing
    % viscosity (SVV) operator applied that prevents any damping below a
    % cutoff wavenumber. The SVV operator adjusts the wavenumbers being
    % damped depending on whether anti-aliasing is applied.
    %
    %
    % - Topic: Initializing
    % - Declaration: WVNonlinearFlux < [WVForcingFluxOperation](/classes/wvnonlinearfluxoperation/)
    properties
        damp

        k_damp % wavenumber at which the significant scale damping starts.
        k_no_damp % wavenumber below which there is zero damping
    end

    methods
        function self = WVAdaptiveViscosity(wvt)
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
            end
            self@WVForcing(wvt,"adaptive svv",WVForcingType(["Spectral","PVSpectral"]));
            self.wvt = wvt;
            self.isClosure = true;
            self.buildDampingOperator();
        end

        function buildDampingOperator(self)
            [K,L,J] = self.wvt.kljGrid;
            M = J*pi/self.wvt.Lz;
            [Qkl,Qj,self.k_no_damp,self.k_damp] = self.spectralVanishingViscosityFilter(shouldAssumeAntialiasing=0);
            prefactor = self.wvt.effectiveHorizontalGridResolution/(pi^2);
            self.damp = -(prefactor*Qkl.*(K.^2 +L.^2));
        end

        function [Qkl,Qj,kl_cutoff, kl_damp] = spectralVanishingViscosityFilter(self,options)
            arguments
                self WVAdaptiveViscosity {mustBeNonempty}
                options.shouldAssumeAntialiasing double {mustBeMember(options.shouldAssumeAntialiasing,[0 1])} = 1
            end
            wvt_ = self.wvt;
            % Builds the spectral vanishing viscosity operator
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
            dampingTimeScale = 1/max(abs(self.damp(:)));
        end
        
        function [Fp, Fm, F0] = addSpectralForcing(self, wvt, Fp, Fm, F0)
            uvMax = wvt.uvMax;
            Fp = Fp + uvMax * self.damp .* wvt.Ap;
            Fm = Fm + uvMax * self.damp .* wvt.Am;
            F0 = F0 + uvMax * self.damp .* wvt.A0;
        end

        function F0 = addPotentialVorticitySpectralForcing(self, wvt, F0)
            F0 = F0 + wvt.uvMax * self.damp .* wvt.A0;
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            force = WVAdaptiveViscosity(wvtX2);
        end
    end
    methods (Static)
        function vars = classRequiredPropertyNames()
            vars = {};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
        end
    end
end