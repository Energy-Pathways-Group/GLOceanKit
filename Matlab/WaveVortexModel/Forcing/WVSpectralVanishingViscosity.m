classdef WVSpectralVanishingViscosity < WVForcing
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
        nu_xy = 0
        nu_z = 0
        damp

        k_damp % wavenumber at which the significant scale damping starts.
        k_no_damp % wavenumber below which there is zero damping
    end
    properties (Dependent)
        uv_damp
    end

    methods
        function self = WVSpectralVanishingViscosity(wvt,options)
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
                options.uv_damp (1,1) double 
                options.w_damp (1,1) double % characteristic speed used to set the damping. Try using wMax
                options.nu_xy (1,1) double
                options.nu_z (1,1) double 
            end
            self@WVForcing(wvt,"spectral vanishing viscosity",WVForcingType(["Spectral","PVSpectral"]));
            self.wvt = wvt;
            self.isClosure = true;

            if isfield(options,'nu_xy')
                self.nu_xy = options.nu_xy;
            else
                if isfield(options,'uv_damp')
                    self.nu_xy = wvt.effectiveHorizontalGridResolution*options.uv_damp/(pi^2);
                else
                    self.nu_xy = 0;
                end
            end

            if isfield(options,'nu_z')
                self.nu_z = options.nu_z;
            else
                if isfield(options,'w_damp')
                    if wvt.shouldAntialias == 1
                        self.nu_z = (3/2)*(wvt.z(2)-wvt.z(1))*options.w_damp/(pi^2);
                    else
                        self.nu_z = (wvt.z(2)-wvt.z(1))*options.w_damp/(pi^2);
                    end
                else
                    self.nu_z = 0;
                end
            end
        end

        function set.nu_xy(self,value)
            self.nu_xy = value;
            self.buildDampingOperator();
        end

        function set.uv_damp(self,uv_damp)
            self.nu_xy = self.wvt.effectiveHorizontalGridResolution*uv_damp/(pi^2);
        end

        function val = get.uv_damp(self)
            val = self.nu_xy*(pi^2)/self.wvt.effectiveHorizontalGridResolution;
        end

        function buildDampingOperator(self)
            [K,L,J] = self.wvt.kljGrid;
            M = J*pi/self.wvt.Lz;
            [Qkl,Qj,self.k_no_damp,self.k_damp] = self.spectralVanishingViscosityFilter(shouldAssumeAntialiasing=0);
            self.damp = -(self.nu_z*Qj.*M.^2 + self.nu_xy*Qkl.*(K.^2 +L.^2));
        end

        function [Qkl,Qj,kl_cutoff, kl_damp] = spectralVanishingViscosityFilter(self,options)
            arguments
                self WVSpectralVanishingViscosity {mustBeNonempty}
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
            Fp = Fp + self.damp .* wvt.Ap;
            Fm = Fm + self.damp .* wvt.Am;
            F0 = F0 + self.damp .* wvt.A0;
        end

        function F0 = addPotentialVorticitySpectralForcing(self, wvt, F0)
            F0 = F0 + self.damp .* wvt.A0;
        end

        function force = forcingWithResolutionOfTransform(self,wvtX2)
            ratio_h = self.wvt.effectiveHorizontalGridResolution/wvtX2.effectiveHorizontalGridResolution;
            ratio_v = self.wvt.effectiveVerticalGridResolution/wvtX2.effectiveVerticalGridResolution;
            force = WVSpectralVanishingViscosity(wvtX2,nu_xy=self.nu_xy/ratio_h,nu_z=self.nu_z/ratio_v);
        end
    end

    methods (Static)
        function vars = classRequiredPropertyNames()
            vars = {'nu_xy','nu_z'};
        end

        function propertyAnnotations = classDefinedPropertyAnnotations()
            arguments (Output)
                propertyAnnotations CAPropertyAnnotation
            end
            propertyAnnotations = CAPropertyAnnotation.empty(0,0);
            propertyAnnotations(end+1) = CANumericProperty('nu_xy', {}, 'm^2 s^{-1}','horizontal viscosity');
            propertyAnnotations(end+1) = CANumericProperty('nu_z', {}, 'm^2 s^{-1}','vertical viscosity');
        end
    end

end