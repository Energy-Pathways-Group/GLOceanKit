classdef SingleModeQGPVE < NonlinearFluxOperation
    properties
        PVA0 % conversion from A0 to PV
        A0PV % conversion from PV to A0
        beta
        damp
        nu
        r
    end

    methods
        function self = SingleModeQGPVE(wvt,options)
            arguments
                wvt WaveVortexTransform {mustBeNonempty}
                options.shouldUseBeta double {mustBeMember(options.shouldUseBeta,[0 1])} = 0 
                options.u_damp (1,1) double = 0.25 % characteristic speed used to set the damping. Try using u_max.
                options.r (1,1) double = 0
                options.fluxName char = 'SingleModeQGPVE'
                options.nu (1,1) double
                options.stateVariables StateVariable = StateVariable.empty()
            end
            fluxVar(1) = StateVariable('F0',{'k','l','j'},'m/s', 'non-linear flux into A0');
            fluxVar(2) = StateVariable('u',{'x','y','z'},'m/s', 'geostrophic velocity x-direction');
            fluxVar(3) = StateVariable('v',{'x','y','z'},'m/s', 'geostrophic velocity y-direction');
%             fluxVar(4) = StateVariable('qgpv',{'x','y','z'},'m/s', 'quasigeostrophic potential vorticity');
            fluxVar = cat(2,fluxVar,options.stateVariables);

            self@NonlinearFluxOperation(options.fluxName,fluxVar);
            self.doesFluxAp = 0;
            self.doesFluxAm = 0;
            self.doesFluxA0 = 1;
            
            AA = ~(wvt.MaskForAliasedModes(jFraction=1));
            self.PVA0 = - wvt.Omega .* wvt.Omega / (wvt.h * wvt.f0);
            self.A0PV = AA./self.PVA0;
            
            % Components to the damping operator (which will multiply A0):
            % 1. Convert A0 to a velocity streamfunction (g/f0)
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
            [K,L] = ndgrid(wvt.k,wvt.l,wvt.j);
            self.damp = -self.r*(K.^2 +L.^2) - self.nu*Qkl.*(-(K.^2 +L.^2)).^2; 
            self.damp = (wvt.g*wvt.h ./(wvt.Omega .* wvt.Omega)) .* AA.* self.damp; % (g/f0) converts A0 into a velocity
            if options.shouldUseBeta == 1
                self.beta = 2 * 7.2921E-5 * cos( wvt.latitude*pi/180. ) / 6.371e6;
            else
                self.beta = 0;
            end
        end

        function dampingTimeScale = dampingTimeScale(self)
            dampingTimeScale = 1/max(max(abs(self.damp)));
        end

        function varargout = Compute(self,wvt,varargin)
            % Apply operator S---defined in (C4) in the manuscript
            Ubar = wvt.UA0 .* wvt.A0;
            Vbar = wvt.VA0 .* wvt.A0;
            PVbar = self.PVA0 .* wvt.A0;

            u_g = wvt.transformToSpatialDomainWithF(Ubar);
            v_g = wvt.transformToSpatialDomainWithF(Vbar);
%             QGPV = wvt.transformToSpatialDomainWithF(PVbar);
            PVx = wvt.transformToSpatialDomainWithF(sqrt(-1)*wvt.k.*PVbar);
            PVy = wvt.transformToSpatialDomainWithF(sqrt(-1)*shiftdim(wvt.l,-1).*PVbar);

            PVnl = u_g.*PVx + v_g.*(PVy+self.beta);
            F0 = -self.A0PV .* wvt.transformFromSpatialDomainWithF(PVnl) + self.damp .* wvt.A0;
            varargout = {F0,u_g,v_g};
        end

        function writeToFile(self,ncfile,wvt)
            arguments
                self NonlinearFluxOperation {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                wvt WaveVortexTransform {mustBeNonempty}
            end
            ncfile.addAttribute('beta',self.beta)
            ncfile.addAttribute('r',self.r)
            ncfile.addAttribute('nu',self.nu)

            attributes = containers.Map();
            attributes('units') = '1/s';
            attributes('description') = 'Linear damping operator applied to A0 to produce damping flux.';
            ncfile.initComplexVariable('L_damp',{'k','l','j'},attributes,'NC_DOUBLE');
            ncfile.setVariable('L_damp',self.damp);
        end

        function nlFlux = nonlinearFluxWithDoubleResolution(self,wvtX2)
            nlFlux = SingleModeQGPVE(wvtX2,r=self.r,shouldUseBeta=(self.beta>0),nu=self.nu/2);
        end

        function flag = isequal(self,other)
            arguments
                self NonlinearFluxOperation
                other NonlinearFluxOperation
            end
            flag = isequal@NonlinearFluxOperation(self,other);
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
                wvt WaveVortexTransform {mustBeNonempty}
            end
            nlFlux = SingleModeQGPVE(wvt,r=ncfile.attributes('r'),nu=ncfile.attributes('nu'),shouldUseBeta=(ncfile.attributes('beta')>0) );
        end
    end

end