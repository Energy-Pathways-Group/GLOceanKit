classdef WVNonlinearFluxQGForced < WVNonlinearFluxQG
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
        EMA0 % energy mask, A0
        FA0 = []  % forcing value, A0
        FTA0  % forcing time, A0 (scalar)
    end

    methods
        function self = WVNonlinearFluxQGForced(wvt,options,newOptions)
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
                options.fluxName char = 'WVNonlinearFluxQGForced'
                options.nu (1,1) double
                options.stateVariables WVVariableAnnotation = WVVariableAnnotation.empty()

                newOptions.EMA0
                newOptions.FA0 = []
                newOptions.FTA0 = 0
            end

            qgArgs = namedargs2cell(options);
            self@WVNonlinearFluxQG(wvt,qgArgs{:});

            self.setGeostrophicForcingCoefficients(zeros(wvt.Nk,wvt.Nl,wvt.Nj),ones(wvt.Nk,wvt.Nl,wvt.Nj),newOptions.FTA0);
        end

        function setGeostrophicForcingCoefficients(self,forcedCoefficients,coefficientMask,relaxationTime)
            arguments
                self WVNonlinearFluxQGForced {mustBeNonempty}
                forcedCoefficients (:,:,:) double {mustBeNonempty}
                coefficientMask (:,:,:) double {mustBeNonempty,mustBeReal}
                relaxationTime (1,1) double = 0
            end

            if relaxationTime > 0
                self.EMA0 = coefficientMask;
                self.FA0 = forcedCoefficients;
                self.FTA0 = relaxationTime;
            else
                % if the relaxation time is zero, then we just want to fix
                % the energy at this mode. So we first make sure the wvt is
                % set to force with those values, then invert the mask, so
                % that it is zeros at the forcing modes.
                self.wvt.A0(logical(coefficientMask)) = forcedCoefficients(logical(coefficientMask));
                self.EMA0 = ~WVTransform.makeHermitian(coefficientMask);
            end
        end

        function varargout = compute(self,wvt,varargin)
            varargout = cell(1,self.nVarOut);
            [varargout{:}] = compute@WVNonlinearFluxQG(self,wvt,varargin{:});

            if self.FTA0 > 0
                varargout{1} = self.EMA0.*(self.FA0 - wvt.A0)/self.FTA0 + varargout{1};
            else
                varargout{1} = self.EMA0 .* varargout{1};
            end
        end

        function writeToFile(self,ncfile,wvt)
            arguments
                self WVNonlinearFluxOperation {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            writeToFile@WVNonlinearFluxQG(self,ncfile,wvt);
            ncfile.addVariable('EMA0',int8(self.EMA0),{'k','l','j'});
            ncfile.addVariable('FA0',self.FA0,{'k','l','j'});
            ncfile.addVariable('FTA0',wvt.FTA0,{});
        end

        function nlFlux = nonlinearFluxWithDoubleResolution(self,wvtX2)
            nlFlux = WVNonlinearFluxQGForced(wvtX2,r=self.r,shouldUseBeta=(self.beta>0),nu=self.nu/2);
        end

        function flag = isequal(self,other)
            arguments
                self WVNonlinearFluxOperation
                other WVNonlinearFluxOperation
            end
            flag = isequal@WVNonlinearFluxQG(self,other);
            flag = flag & isequal(self.EMA0, other.EMA0);
            flag = flag & isequal(self.FA0, other.FA0);
            flag = flag & isequal(self.FTA0, other.FTA0);
        end

    end

    methods (Static)
        function nlFlux = nonlinearFluxFromFile(ncfile,wvt)
            arguments
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            nlFlux = WVNonlinearFluxQGForced(wvt,r=ncfile.attributes('r'),nu=ncfile.attributes('nu'),shouldUseBeta=(ncfile.attributes('beta')>0) );
            nlFlux.EMA0 = logical(ncfile.readVariables('EMA0'));
            nlFlux.FA0 = ncfile.readVariables('FA0');
            nlFlux.FTA0 = ncfile.readVariables('FTA0');
        end
    end

end