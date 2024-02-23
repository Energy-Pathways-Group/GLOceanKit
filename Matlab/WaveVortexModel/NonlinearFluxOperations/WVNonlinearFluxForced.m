classdef WVNonlinearFluxForced < WVNonlinearFlux
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
    % model = WVModel(wvt,nonlinearFlux=WVNonlinearFluxForced(wvt,uv_damp=wvt.uMax));
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
        function self = WVNonlinearFluxForced(wvt,options)
            % initialize WVNonlinearFluxForced
            %
            % - Topic: Initialization
            % - Declaration: nlFlux = WVNonlinearFluxForced(wvt,options)
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
            
            qgArgs = namedargs2cell(options);
            self@WVNonlinearFlux(wvt,qgArgs{:});
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
                self WVNonlinearFluxForced {mustBeNonempty}
                Apbar (:,:,:) double {mustBeNonempty}
                Ambar (:,:,:) double {mustBeNonempty}
                options.MAp (:,:,:) logical = abs(Apbar) > 0
                options.MAm (:,:,:) logical = abs(Ambar) > 0
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

            % multiply by the anti-alias filter so we don't force in the
            % aliased region.
            self.Apbar = self.AA .* Apbar;
            self.MAp = options.MAp;
            self.tauP = options.tauP;

            self.Ambar = self.AA .* Ambar;
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
                self WVNonlinearFluxForced {mustBeNonempty}
                A0bar (:,:,:) double {mustBeNonempty}
                options.MA0 (:,:,:) logical = abs(A0bar) > 0
                options.tau0 (1,1) double = 0
            end

            self.A0bar = self.AA .* A0bar;
            self.MA0 = options.MA0;
            self.tau0 = options.tau0;
        end

        function varargout = compute(self,wvt,varargin)
            varargout = cell(1,self.nVarOut);
            [varargout{:}] = compute@WVNonlinearFlux(self,wvt,varargin{:});

            if self.tauP > 0
                varargout{1} = self.MAp.*(self.Apbar - wvt.Ap)/self.tauP + varargout{1};
            elseif ~isempty(self.MAp)
                varargout{1} = (~self.MAp) .* varargout{1};
            end
            if self.tauM > 0
                varargout{2} = self.MAm.*(self.Ambar - wvt.Am)/self.tauM + varargout{2};
            elseif ~isempty(self.MAm)
                varargout{2} = (~self.MAm) .* varargout{2};
            end
            if self.tau0 > 0
                varargout{3} = self.MA0.*(self.A0bar - wvt.A0)/self.tau0 + varargout{3};
            elseif ~isempty(self.MA0)
                varargout{3} = (~self.MA0) .* varargout{3};
            end
        end

        function writeToFile(self,ncfile,wvt)
            arguments
                self WVNonlinearFluxOperation {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            writeToFile@WVNonlinearFlux(self,ncfile,wvt);
            self.addVariableOfType(ncfile,wvt,'MAp',self.MAp,@(v) uint8(v));
            self.addVariableOfType(ncfile,wvt,'Apbar',self.Apbar,@(v) v);
            ncfile.addVariable('tauP',self.tauP,{});
            self.addVariableOfType(ncfile,wvt,'MAm',self.MAm,@(v) uint8(v));
            self.addVariableOfType(ncfile,wvt,'Ambar',self.Ambar,@(v) v);
            ncfile.addVariable('tauM',self.tauM,{});
            self.addVariableOfType(ncfile,wvt,'MA0',self.MA0,@(v) uint8(v));
            self.addVariableOfType(ncfile,wvt,'A0bar',self.A0bar,@(v) v);
            ncfile.addVariable('tau0',self.tau0,{});
        end

        function addVariableOfType(self,ncfile,wvt,name,var,type)
            if isempty(var)
                ncfile.addVariable(name,type(zeros(wvt.Nk,wvt.Nl,wvt.Nj)),{'k','l','j'});
            else
                ncfile.addVariable(name,type(var),{'k','l','j'});
            end
        end

        function nlFlux = nonlinearFluxWithResolutionOfTransform(self,wvtX2)
            nlFlux = WVNonlinearFluxForced(wvtX2,nu_xy=self.nu_xy/2,nu_z=self.nu_z/2,shouldAntialias=self.shouldAntialias);
            if ~isempty(self.MAp)
                nlFlux.MAp = WVTransform.spectralVariableWithResolution(self.MAp,[wvtX2.Nk wvtX2.Nl wvtX2.Nj]);
            end
            if ~isempty(self.Apbar)
                nlFlux.Apbar = WVTransform.spectralVariableWithResolution(self.Apbar,[wvtX2.Nk wvtX2.Nl wvtX2.Nj]);
            end
            nlFlux.tauP = self.tauP;

            if ~isempty(self.MAm)
                nlFlux.MAm = WVTransform.spectralVariableWithResolution(self.MAm,[wvtX2.Nk wvtX2.Nl wvtX2.Nj]);
            end
            if ~isempty(self.Ambar)
                nlFlux.Ambar = WVTransform.spectralVariableWithResolution(self.Ambar,[wvtX2.Nk wvtX2.Nl wvtX2.Nj]);
            end
            nlFlux.tauM = self.tauM;

            if ~isempty(self.MA0)
                nlFlux.MA0 = WVTransform.spectralVariableWithResolution(self.MA0,[wvtX2.Nk wvtX2.Nl wvtX2.Nj]);
            end
            if ~isempty(self.A0bar)
                nlFlux.A0bar = WVTransform.spectralVariableWithResolution(self.A0bar,[wvtX2.Nk wvtX2.Nl wvtX2.Nj]);
            end
            nlFlux.tau0 = self.tau0;
        end

        function flag = isequal(self,other)
            arguments
                self WVNonlinearFluxOperation
                other WVNonlinearFluxOperation
            end
            flag = isequal@WVNonlinearFlux(self,other);
            flag = flag & isequal(self.MAp, other.MAp);
            flag = flag & isequal(self.Apbar, other.Apbar);
            flag = flag & isequal(self.tauP, other.tauP);
            flag = flag & isequal(self.MAm, other.MAm);
            flag = flag & isequal(self.Ambar, other.Ambar);
            flag = flag & isequal(self.tauM, other.tauM);
            flag = flag & isequal(self.MA0, other.MA0);
            flag = flag & isequal(self.A0bar, other.A0bar);
            flag = flag & isequal(self.tau0, other.tau0);
        end
    end

    methods (Static)
        function nlFlux = nonlinearFluxFromFile(ncfile,wvt)
            arguments
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            nlFlux = WVNonlinearFluxForced(wvt,nu_xy=ncfile.attributes('nu_xy'),nu_z=ncfile.attributes('nu_z'),shouldAntialias=ncfile.attributes('shouldAntialias') );

            nlFlux.MAp = logical(ncfile.readVariables('MAp'));
            nlFlux.Apbar = ncfile.readVariables('Apbar');
            nlFlux.tauP = ncfile.readVariables('tauP');
            if ~any(nlFlux.MAp)
                nlFlux.MAp = [];
                nlFlux.Apbar = [];
            end       
            
            nlFlux.MAm = logical(ncfile.readVariables('MAm'));
            nlFlux.Ambar = ncfile.readVariables('Ambar');
            nlFlux.tauM = ncfile.readVariables('tauM');
            if ~any(nlFlux.MAm)
                nlFlux.MAm = [];
                nlFlux.Ambar = [];
            end

            nlFlux.MA0 = logical(ncfile.readVariables('MA0'));
            nlFlux.A0bar = ncfile.readVariables('A0bar');
            nlFlux.tau0 = ncfile.readVariables('tau0');
            if ~any(nlFlux.MA0)
                nlFlux.MA0 = [];
                nlFlux.A0bar = [];
            end
        end
    end

end