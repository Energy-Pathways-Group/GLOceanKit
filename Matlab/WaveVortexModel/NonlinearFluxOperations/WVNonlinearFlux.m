classdef WVNonlinearFlux < WVNonlinearFluxOperation
    % 3D nonlinear flux for Boussinesq flow, appropriate for numerical modeling
    %
    % Computes the nonlinear flux for a Boussinesq model, and has options
    % for anti-aliasing and damping appropriate for running a numerical
    % model. This is *not* the simplest implementation, but instead adds
    % some complexity in favor of speed. The [BoussinesqSpatial](/classes/boussinesqspatial/) class
    % shows a simple implementation.
    %
    % The damping is a simple Laplacian, but with a spectral vanishing
    % viscosity (SVV) operator applied that prevents any damping below a
    % cutoff wavenumber. The SVV operator adjusts the wavenumbers being
    % damped depending on whether anti-aliasing is applied.
    %
    % This is most often used when initializing a model, e.g.,
    %
    % ```matlab
    % model = WVModel(wvt,nonlinearFlux=WVNonlinearFlux(wvt,shouldAntialias=1,uv_damp=wvt.uMax));
    % ```
    %
    % - Topic: Initializing
    % - Declaration: WVNonlinearFlux < [WVNonlinearFluxOperation](/classes/wvnonlinearfluxoperation/)
    properties
        wvt
        shouldAntialias = 0
        AA
        nu_xy = 0
        nu_z = 0
        r
        damp
        dLnN2 = 0
        beta
        betaA0
        k_damp % wavenumber at which the first small scale damping starts.
    end
    properties (Dependent)
        uv_damp
    end

    methods
        function self = WVNonlinearFlux(wvt,options)
            % initialize the WVNonlinearFlux nonlinear flux
            %
            % - Declaration: nlFlux = WVNonlinearFlux(wvt,options)
            % - Parameter wvt: a WVTransform instance
            % - Parameter uv_damp: (optional) characteristic speed used to set the damping. Try using wvt.uMax.
            % - Parameter w_damp: (optional) characteristic speed used to set the damping. Try using wvt.wMax.
            % - Parameter nu_xy: (optional) coefficient for damping
            % - Parameter nu_z: (optional) coefficient for damping
            % - Parameter shouldAntialias: (optional) a Boolean indicating whether or not to antialias (default 1)
            % - Returns nlFlux: a WVNonlinearFlux instance
            arguments
                wvt WVTransform {mustBeNonempty}
                options.uv_damp (1,1) double 
                options.w_damp (1,1) double % characteristic speed used to set the damping. Try using wMax
                options.nu_xy (1,1) double
                options.nu_z (1,1) double 
                options.r (1,1) double {mustBeNonnegative} = 0 % linear bottom friction, try 1/(200*86400) https://www.nemo-ocean.eu/doc/node70.html
                options.shouldAntialias double = 1
                options.shouldUseBeta double {mustBeMember(options.shouldUseBeta,[0 1])} = 0
            end
            fluxVar(1) = WVVariableAnnotation('Fp',{'k','l','j'},'m/s2', 'non-linear flux into Ap');
            fluxVar(2) = WVVariableAnnotation('Fm',{'k','l','j'},'m/s2', 'non-linear flux into Am');
            fluxVar(3) = WVVariableAnnotation('F0',{'k','l','j'},'m/s', 'non-linear flux into A0');

            self@WVNonlinearFluxOperation('WVNonlinearFlux',fluxVar);
            self.wvt = wvt;
            self.shouldAntialias = options.shouldAntialias;
            self.r = options.r;
            
            if self.shouldAntialias == 1
                self.AA = ~(wvt.maskForAliasedModes(jFraction=2/3));
                wvt.Ap = self.AA .* wvt.Ap;
                wvt.Am = self.AA .* wvt.Am;
                wvt.A0 = self.AA .* wvt.A0;
            else
                self.AA = 1;
            end
            
            if isa(wvt,'WVTransformConstantStratification')
                self.dLnN2 = 0;
            elseif isa(wvt,'WVTransformHydrostatic')
                self.dLnN2 = shiftdim(wvt.dLnN2,-2);
            else
                self.dLnN2 = shiftdim(wvt.dLnN2,-2);
                warning('WVTransform not recognized.')
            end

            if isfield(options,'nu_xy')
                self.nu_xy = options.nu_xy;
            else
                if isfield(options,'uv_damp')
                    if self.shouldAntialias == 1
                        self.nu_xy = (3/2)*(wvt.x(2)-wvt.x(1))*options.uv_damp/(pi^2);
                    else
                        self.nu_xy = (wvt.x(2)-wvt.x(1))*options.uv_damp/(pi^2);
                    end
                else
                    self.nu_xy = 0;
                end
            end

            if isfield(options,'nu_z')
                self.nu_z = options.nu_z;
            else
                if isfield(options,'w_damp')
                    if self.shouldAntialias == 1
                        self.nu_z = (3/2)*(wvt.z(2)-wvt.z(1))*options.w_damp/(pi^2);
                    else
                        self.nu_z = (wvt.z(2)-wvt.z(1))*options.w_damp/(pi^2);
                    end
                else
                    self.nu_z = 0;
                end
            end
            
            if options.shouldUseBeta == 1
                self.beta = 2 * 7.2921E-5 * cos( wvt.latitude*pi/180. ) / 6.371e6;
                self.betaA0 = -self.beta * (wvt.VA0 ./ wvt.A0_QGPV_factor);
                self.betaA0(1,1,1) = 0;
            else
                self.beta = 0;
                self.betaA0 = 0;
            end

            [K,L,J] = ndgrid(wvt.k,wvt.l,wvt.j);
            M = J*pi/wvt.Lz;
            self.damp = -(self.nu_z*M.^2 + self.nu_xy*(K.^2 +L.^2));

            Qkl = wvt.spectralVanishingViscosityFilter(shouldAssumeAntialiasing=self.shouldAntialias);
            self.damp = Qkl.*self.damp;
        end

        function set.nu_xy(self,value)
            self.nu_xy = value;
            self.buildDampingOperator();
        end

        function set.uv_damp(self,uv_damp)
            if self.shouldAntialias == 1
                self.nu_xy = (3/2)*(self.wvt.x(2)-self.wvt.x(1))*uv_damp/(pi^2);
            else
                self.nu_xy = (self.wvt.x(2)-self.wvt.x(1))*uv_damp/(pi^2);
            end
        end

        function val = get.uv_damp(self)
            if self.shouldAntialias == 1
                val = self.nu_xy*(pi^2)*(2/3)/(self.wvt.x(2)-self.wvt.x(1));
            else
                val = self.nu_xy*(pi^2)/(self.wvt.x(2)-self.wvt.x(1));
            end
        end

        function buildDampingOperator(self)
            [K,L,J] = self.wvt.kljGrid;
            M = J*pi/self.wvt.Lz;
            self.damp = -(self.nu_z*M.^2 + self.nu_xy*(K.^2 +L.^2));
            [Qkl,~,self.k_damp] = self.wvt.spectralVanishingViscosityFilter(shouldAssumeAntialiasing=self.shouldAntialias);
            self.damp = Qkl.*self.damp;
        end

        function dampingTimeScale = dampingTimeScale(self)
            dampingTimeScale = 1/max(abs(self.damp(:)));
        end
        
        function varargout = spatialFlux(self,wvt,varargin)
            % a subclass can override this, and then modify the spatial
            % fluxes that get returned.
            phase = exp(wvt.iOmega*(wvt.t-wvt.t0));
            Apt = wvt.Ap .* phase;
            Amt = wvt.Am .* conj(phase);
            A0t = wvt.A0;

            % Apply operator S---defined in (C4) in the manuscript
            Ubar = wvt.UAp.*Apt + wvt.UAm.*Amt + wvt.UA0.*A0t;
            Vbar = wvt.VAp.*Apt + wvt.VAm.*Amt + wvt.VA0.*A0t;
            Wbar = wvt.WAp.*Apt + wvt.WAm.*Amt;
            Nbar = wvt.NAp.*Apt + wvt.NAm.*Amt + wvt.NA0.*A0t;

            % Finishing applying S, but also compute derivatives at the
            % same time
            [U,Ux,Uy,Uz] = wvt.transformToSpatialDomainWithFAllDerivatives(Ubar);
            [V,Vx,Vy,Vz] = wvt.transformToSpatialDomainWithFAllDerivatives(Vbar);
            W = wvt.transformToSpatialDomainWithG(Wbar);
            [ETA,ETAx,ETAy,ETAz] = wvt.transformToSpatialDomainWithGAllDerivatives(Nbar);

            % compute the nonlinear terms in the spatial domain
            % (pseudospectral!)
            uNL = -U.*Ux - V.*Uy - W.*Uz;
            vNL = -U.*Vx - V.*Vy - W.*Vz;
            nNL = -U.*ETAx - V.*ETAy - W.*(ETAz + ETA.*self.dLnN2);
            
            % bottom friction
            uNL(:,:,1) = uNL(:,:,1) - self.r*U(:,:,1);
            vNL(:,:,1) = vNL(:,:,1) - self.r*V(:,:,1);
            
            if nargout == 5
                uvw = struct('u',U,'v',V,'w',W,'eta',ETA,'ux',Ux,'uy',Uy,'uz',Uz,'vx',Vx,'vy',Vy,'vz',Vz,'etax',ETAx,'etay',ETAy,'etaz',ETAz);
                varargout = {uNL,vNL,nNL,phase,uvw};
            else
                varargout = {uNL,vNL,nNL,phase};
            end
        end

        function varargout = compute(self,wvt,varargin)
            [uNL,vNL,nNL,phase] = self.spatialFlux(wvt,varargin);

            % Now apply the operator S^{-1} and then T_\omega^{-1}
            uNLbar = wvt.transformFromSpatialDomainWithF(uNL);
            vNLbar = wvt.transformFromSpatialDomainWithF(vNL);
            nNLbar = wvt.transformFromSpatialDomainWithG(nNL);

            Fp = self.AA .* (self.damp .* wvt.Ap + (wvt.ApU.*uNLbar + wvt.ApV.*vNLbar + wvt.ApN.*nNLbar) .* conj(phase));
            Fm = self.AA .* (self.damp .* wvt.Am + (wvt.AmU.*uNLbar + wvt.AmV.*vNLbar + wvt.AmN.*nNLbar) .* phase);
            F0 = self.AA .* ((self.damp + self.betaA0) .* wvt.A0 + (wvt.A0U.*uNLbar + wvt.A0V.*vNLbar + wvt.A0N.*nNLbar));

            varargout = {Fp,Fm,F0};
        end

        function writeToFile(self,ncfile,wvt)
            arguments
                self WVNonlinearFluxOperation {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            ncfile.addAttribute('nu_xy',self.nu_xy)
            ncfile.addAttribute('nu_z',self.nu_z)
            ncfile.addAttribute('r',self.r)
            ncfile.addAttribute('shouldAntialias',self.shouldAntialias)
        end

        function nlFlux = nonlinearFluxWithResolutionOfTransform(self,wvtX2)
            nlFlux = WVNonlinearFlux(wvtX2,nu_xy=self.nu_xy/2,nu_z=self.nu_z/2,shouldAntialias=self.shouldAntialias);
        end

        function flag = isequal(self,other)
            arguments
                self WVNonlinearFluxOperation
                other WVNonlinearFluxOperation
            end
            flag = isequal@WVNonlinearFluxOperation(self,other);
            flag = flag & isequal(self.shouldAntialias, other.shouldAntialias);
            flag = flag & isequal(self.AA, other.AA);
            flag = flag & isequal(self.nu_xy, other.nu_xy);
            flag = flag & isequal(self.nu_z,other.nu_z);
            flag = flag & isequal(self.r,other.r);
            flag = flag & isequal(self.damp, other.damp);
        end
    end

    methods (Static)
        function nlFlux = nonlinearFluxFromFile(ncfile,wvt)
            arguments
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            nlFlux = WVNonlinearFlux(wvt,nu_xy=ncfile.attributes('nu_xy'),nu_z=ncfile.attributes('nu_z'),r=ncfile.attributes('r'),shouldAntialias=ncfile.attributes('shouldAntialias') );
        end
    end

end