classdef WVSpectralVanishingViscosity < handle
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
        wvt
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

            [K,L,J] = self.wvt.kljGrid;
            M = J*pi/wvt.Lz;
            self.damp = -(self.nu_z*M.^2 + self.nu_xy*(K.^2 +L.^2));

            Qkl = wvt.spectralVanishingViscosityFilter(shouldAssumeAntialiasing=0);
            self.damp = Qkl.*self.damp;
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
            [Qkl,Qj,self.k_no_damp,self.k_damp] = self.wvt.spectralVanishingViscosityFilter(shouldAssumeAntialiasing=0);
            self.damp = -(self.nu_z*Qj.*M.^2 + self.nu_xy*Qkl.*(K.^2 +L.^2));
        end

        function dampingTimeScale = dampingTimeScale(self)
            dampingTimeScale = 1/max(abs(self.damp(:)));
        end
        
        function [Fp, Fm, F0] = addForcing(self, wvt, Fp, Fm, F0)
            Fp = Fp + self.damp .* wvt.Ap;
            Fm = Fm + self.damp .* wvt.Am;
            F0 = F0 + self.damp .* wvt.A0;
        end

        function writeToFile(self,group,wvt)
            arguments
                self WVForcingFluxOperation {mustBeNonempty}
                group NetCDFGroup {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            group.addAttribute('nu_xy',self.nu_xy)
            group.addAttribute('nu_z',self.nu_z)
        end

        function forcingFlux = forcingFluxWithResolutionOfTransform(self,wvtX2)
            ratio_h = self.wvt.effectiveHorizontalGridResolution/wvtX2.effectiveHorizontalGridResolution;
            ratio_v = self.wvt.effectiveVerticalGridResolution/wvtX2.effectiveVerticalGridResolution;
            forcingFlux = WVSpectralVanishingViscosity(wvtX2,nu_xy=self.nu_xy/ratio_h,nu_z=self.nu_z/ratio_v);
        end

        function flag = isequal(self,other)
            arguments
                self WVForcingFluxOperation
                other WVForcingFluxOperation
            end
            flag = isequal@WVForcingFluxOperation(self,other);
            flag = flag & isequal(self.nu_xy, other.nu_xy);
            flag = flag & isequal(self.nu_z,other.nu_z);
            flag = flag & isequal(self.damp, other.damp);
        end
    end

    methods (Static)
        function forcingFlux = forcingFluxFromFile(ncfile,wvt)
            arguments
                ncfile NetCDFGroup {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            forcingFlux = WVSpectralVanishingViscosity(wvt,nu_xy=ncfile.attributes('nu_xy'),nu_z=ncfile.attributes('nu_z') );
        end
    end

end