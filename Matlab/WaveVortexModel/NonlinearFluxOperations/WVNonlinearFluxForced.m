classdef WVNonlinearFluxForced < WVNonlinearFlux
    % 3D forced nonlinear flux for Boussinesq flow
    %
    properties
        MA0         % Forcing mask, A0. 1s at the forced modes, 0s at the unforced modes
        MAp         % Forcing mask, Ap. 1s at the forced modes, 0s at the unforced modes
        MAm         % Forcing mask, Am. 1s at the forced modes, 0s at the unforced modes

        A0bar = []  % A0 'mean' value to relax to
        Apbar = []  % Ap 'mean' value to relax to
        Ambar = []  % Am 'mean' value to relax to

        tau0        % A0 relaxation time
        tauP        % Ap relaxation time
        tauM        % Am relaxation time
    end

    methods
        function self = WVNonlinearFluxForced(wvt,options)
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
                options.shouldAntialias double = 1
            end
            
            qgArgs = namedargs2cell(options);
            self@WVNonlinearFlux(wvt,qgArgs{:});
        end

        function varargout = compute(self,wvt,varargin)
            varargout = cell(1,self.nVarOut);
            [varargout{:}] = compute@WVNonlinearFlux(self,wvt,varargin{:});

            if self.tau0P > 0
                varargout{1} = self.MAp.*(self.Apbar - wvt.Ap)/self.tau0P + varargout{1};
            else
                varargout{1} = (~self.MAp) .* varargout{1};
            end
            if self.tau0M > 0
                varargout{2} = self.MAm.*(self.Ambar - wvt.Am)/self.tau0M + varargout{2};
            else
                varargout{2} = (~self.MAm) .* varargout{2};
            end
            if self.tau0 > 0
                varargout{3} = self.MA0.*(self.A0bar - wvt.A0)/self.tau0 + varargout{3};
            else
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
            ncfile.addVariable('MAp',int8(self.MAp),{'k','l','j'});
            ncfile.addVariable('Apbar',self.Apbar,{'k','l','j'});
            ncfile.addVariable('tauP',wvt.tauP,{});
            ncfile.addVariable('MAm',int8(self.MAm),{'k','l','j'});
            ncfile.addVariable('Ambar',self.Ambar,{'k','l','j'});
            ncfile.addVariable('tauM',wvt.tauM,{});
            ncfile.addVariable('MA0',int8(self.MA0),{'k','l','j'});
            ncfile.addVariable('A0bar',self.A0bar,{'k','l','j'});
            ncfile.addVariable('tau0',wvt.tau0,{});
        end

        function nlFlux = nonlinearFluxWithDoubleResolution(self,wvtX2)
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
            nlFlux.MAm = logical(ncfile.readVariables('MAm'));
            nlFlux.Ambar = ncfile.readVariables('Ambar');
            nlFlux.tauM = ncfile.readVariables('tauM');
            nlFlux.MA0 = logical(ncfile.readVariables('MA0'));
            nlFlux.A0bar = ncfile.readVariables('A0bar');
            nlFlux.tau0 = ncfile.readVariables('tau0');
        end
    end

end