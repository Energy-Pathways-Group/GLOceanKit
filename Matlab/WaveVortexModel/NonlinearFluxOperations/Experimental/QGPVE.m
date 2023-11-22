classdef QGPVE < WVNonlinearFluxOperation
    % 3D quasigeostrophic potential vorticity flux
    %
    % The 3D quasigeostrophic potential vorticity flux will only use and
    % modify the A0 coefficients.
    %
    % To initialize the QGPVE,
    %
    % ```matlab
    % model = WVModel(wvt,nonlinearFlux=QGPVE(wvt,shouldUseBeta=1,u_damp=wvt.uMax));
    % ```
    %
    % - Topic: Initializing
    % - Declaration: QGPVE < [WVNonlinearFluxOperation](/classes/wvnonlinearfluxoperation/)
    properties
        PVA0 % conversion from A0 to PV
        A0PV % conversion from PV to A0
        RVA0 % conversion from A0 to RV
        beta
        damp
        nu
        r
    end

    methods
        function self = QGPVE(wvt,options)
            % initialize 3D quasigeostrophic potential vorticity flux
            %
            % - Declaration: nlFlux = QGPVE(wvt,options)
            % - Parameter wvt: a WVTransform instance
            % - Parameter shouldUseBeta: (optional) a Boolean indicating whether or not to include beta in the flux
            % - Parameter u_damp: (optional) characteristic speed used to set the damping. Try using wvt.uMax
            % - Parameter r: (optional) bottom friction
            % - Parameter nu: (optional) coefficient for damping
            % - Returns nlFlux: a QGPVE instance
            arguments
                wvt WVTransform {mustBeNonempty}
                options.shouldUseBeta double {mustBeMember(options.shouldUseBeta,[0 1])} = 0 
                options.u_damp (1,1) double = 0.25 % characteristic speed used to set the damping. Try using uMax.
                options.r (1,1) double = 0
                options.fluxName char = 'QGPVE'
                options.nu (1,1) double
                options.stateVariables WVVariableAnnotation = WVVariableAnnotation.empty()
            end
            fluxVar(1) = WVVariableAnnotation('F0',{'k','l','j'},'m/s', 'non-linear flux into A0');
            fluxVar(2) = WVVariableAnnotation('u',{'x','y','z'},'m/s', 'geostrophic velocity x-direction');
            fluxVar(2).attributes('standard_name') = 'eastward_sea_water_velocity';
            fluxVar(3) = WVVariableAnnotation('v',{'x','y','z'},'m/s', 'geostrophic velocity y-direction');
            fluxVar(3).attributes('standard_name') = 'northward_sea_water_velocity';
            fluxVar = cat(2,fluxVar,options.stateVariables);

            self@WVNonlinearFluxOperation(options.fluxName,fluxVar);
            self.doesFluxAp = 0;
            self.doesFluxAm = 0;
            self.doesFluxA0 = 1;
            
            AA = ~(wvt.maskForAliasedModes(jFraction=2/3));
            self.PVA0 = wvt.A0_QGPV_factor;
            self.RVA0 = - wvt.g * (wvt.Kh .* wvt.Kh) / wvt.f;
            self.A0PV = AA./self.PVA0;
            self.A0PV(1,1,1) = 0;
            
            % Components to the damping operator (which will multiply A0):
            % 1. Convert A0 to a velocity streamfunction (g/f)
            % 2. bi-harmonic operator
            % 3. spectral vanish viscosity
            % 4. nu is set to create approximately Reynolds number=1.
            Qkl = wvt.spectralVanishingViscosityFilter(shouldAssumeAntialiasing=1);
            if isfield(options,"nu")
                self.nu = options.nu;
            else
                self.nu = (3/2)*(wvt.x(2)-wvt.x(1))*options.u_damp/(pi^2);
            end
            self.r = options.r;
            [K,L] = wvt.kljGrid;
            if wvt.isBarotropic
                friction = -self.r*(K.^2 +L.^2);
            else
                % only apply friction to the bottom
                friction = -0*self.r*(K.^2 +L.^2);
                mask = zeros(length(wvt.z),1);
                mask(1) = 1;
                maskTransform = abs(wvt.FMatrix * mask);
                friction = friction.*shiftdim(maskTransform,-2);
            end
            self.damp = friction - self.nu*Qkl.*(-(K.^2 +L.^2)).^2; 
            self.damp = -(wvt.g/wvt.f) * self.A0PV .* self.damp; % (g/f) converts A0 into a velocity
            if options.shouldUseBeta == 1
                self.beta = 2 * 7.2921E-5 * cos( wvt.latitude*pi/180. ) / 6.371e6;
            else
                self.beta = 0;
            end
        end

        function dampingTimeScale = dampingTimeScale(self)
            dampingTimeScale = 1/max(abs(self.damp(:)));
        end

        function varargout = compute(self,wvt,varargin)
            % Apply operator S---defined in (C4) in the manuscript
            Ubar = wvt.UA0 .* wvt.A0;
            Vbar = wvt.VA0 .* wvt.A0;
            PVbar = self.PVA0 .* wvt.A0;

            u_g = wvt.transformToSpatialDomainWithF(Ubar);
            v_g = wvt.transformToSpatialDomainWithF(Vbar);
            PVx = wvt.transformToSpatialDomainWithF(sqrt(-1)*wvt.k.*PVbar);
            PVy = wvt.transformToSpatialDomainWithF(sqrt(-1)*shiftdim(wvt.l,-1).*PVbar);

            if wvt.isBarotropic
                PVnl = u_g.*PVx + v_g.*(PVy+self.beta);
            else
                mask = zeros(size(wvt.X));
                mask(:,:,1) = 1;
                PVabs = -self.r * mask.*wvt.transformToSpatialDomainWithF(self.RVA0 .* wvt.A0);
                PVnl = u_g.*PVx + v_g.*(PVy+self.beta) + PVabs;
            end
            F0 = -self.A0PV .* wvt.transformFromSpatialDomainWithF(PVnl) + self.damp .* wvt.A0;
            varargout = {F0,u_g,v_g};
        end

        function writeToFile(self,ncfile,wvt)
            arguments
                self WVNonlinearFluxOperation {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            ncfile.addAttribute('beta',self.beta)
            ncfile.addAttribute('r',self.r)
            ncfile.addAttribute('nu',self.nu)

            attributes = containers.Map();
            attributes('units') = '1/s';
            attributes('long_name') = 'Linear damping operator applied to A0 to produce damping flux.';
            ncfile.initComplexVariable('L_damp',{'k','l','j'},attributes,'NC_DOUBLE');
            ncfile.setVariable('L_damp',self.damp);
        end

        function nlFlux = nonlinearFluxWithResolutionOfTransform(self,wvtX2)
            nlFlux = QGPVE(wvtX2,r=self.r,shouldUseBeta=(self.beta>0),nu=self.nu/2);
        end

        function flag = isequal(self,other)
            arguments
                self WVNonlinearFluxOperation
                other WVNonlinearFluxOperation
            end
            flag = isequal@WVNonlinearFluxOperation(self,other);
            flag = flag & isequal(self.PVA0, other.PVA0);
            flag = flag & isequal(self.beta, other.beta);
            flag = flag & isequal(self.damp, other.damp);
            flag = flag & isequal(self.nu,other.nu);
            flag = flag & isequal(self.r, other.r);
        end

    end

    methods (Static)
        function nlFlux = nonlinearFluxFromFile(ncfile,wvt)
            arguments
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            nlFlux = QGPVE(wvt,r=ncfile.attributes('r'),nu=ncfile.attributes('nu'),shouldUseBeta=(ncfile.attributes('beta')>0) );
        end
    end

end