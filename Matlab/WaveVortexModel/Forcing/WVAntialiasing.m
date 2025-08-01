classdef WVAntialiasing < WVForcing
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
        M
    end

    methods
        function self = WVAntialiasing(wvt,options)
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
                options.Nj
            end
            self@WVForcing(wvt,"antialias filter",WVForcingType(["Spectral","PVSpectral"]));
            self.priority = 127;
            if wvt.shouldAntialias
                error("WVAntialiasing:AntialiasingNotSupported","Antialiasing is not supported for a transform that aliases at the transform level.");
            end
            self.wvt = wvt;
            self.isClosure = true;
            Aklz = WVGeometryDoublyPeriodic.maskForAliasedModes(wvt.k_dft,wvt.l_dft,wvt.Nj);
            self.M = wvt.transformFromDFTGridToWVGrid(Aklz);
            if ~isfield(options,"Nj")
                options.Nj = floor(2*wvt.Nj/3);
            end
            self.M(wvt.J > (options.Nj-1)) = 1;
        end

        function effectiveHorizontalGridResolution = effectiveHorizontalGridResolution(self)
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
                self WVAntialiasing
            end
            effectiveHorizontalGridResolution = pi/max(max(abs(self.wvt.L(~self.M)),abs(self.wvt.K(~self.M))));
        end

        function j_max = effectiveJMax(self)
            arguments
                self WVAntialiasing
            end
            j_max = max(abs(self.wvt.J(~self.M)));
        end
        
        function [Fp, Fm, F0] = addSpectralForcing(self, wvt, Fp, Fm, F0)
            Fp = Fp - self.M .* Fp;
            Fm = Fm - self.M .* Fm;
            F0 = F0 - self.M .* F0;
        end

        function F0 = addPotentialVorticitySpectralForcing(self, wvt, F0)
            F0 = F0 - self.M .* F0;
        end

        function [Ap, Am, A0] = setSpectralAmplitude(self, wvt, Ap, Am, A0)
            Ap = (~self.M) .* Ap;
            Am = (~self.M) .* Am;
            A0 = (~self.M) .* A0;
        end



        function A0 = setPotentialVorticitySpectralAmplitude(self, wvt, A0)
            arguments
                self WVForcing
                wvt WVTransform
                A0 
            end
            A0 = (~self.M) .* A0;
        end

        function F0 = setPotentialVorticitySpectralForcing(self, wvt, F0)
            arguments
                self WVForcing
                wvt WVTransform
                F0
            end
            F0 = F0 - self.M .* F0;
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            force = WVAntialiasing(wvtX2);
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