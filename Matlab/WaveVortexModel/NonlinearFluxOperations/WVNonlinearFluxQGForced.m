classdef WVNonlinearFluxQGForced < WVNonlinearFluxQG
    % 3D forced quasigeostrophic potential vorticity flux
    %
    % The 3D quasigeostrophic potential vorticity flux will only use and
    % modify the A0 coefficients.
    %
    % $$
    % \frac{\partial}{\partial t} A_0^{klj} = \underbrace{M_{A_0}^{klj} \left(\bar{A}_0^{klj}  - A_0^{klj} \right)/ \tau_0}_{F_\textrm{force}} + F_0^{klj} + F_\textrm{damp}^{klj}
    % $$
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
        MA0         % Forcing mask, A0. 1s at the forced modes, 0s at the unforced modes
        A0bar = []  % A0 'mean' value to relax to
        tau0        % relaxation time
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

                newOptions.MA0
                newOptions.A0bar = []
                newOptions.tau0 = 0
            end

            qgArgs = namedargs2cell(options);
            self@WVNonlinearFluxQG(wvt,qgArgs{:});

            % self.setGeostrophicForcingCoefficients(zeros(wvt.Nk,wvt.Nl,wvt.Nj),ones(wvt.Nk,wvt.Nl,wvt.Nj),newOptions.FTA0);
        end

        function setGeostrophicForcingCoefficients(self,A0bar,options)
            % set forcing values for the geostrophic part of the flow
            %
            % $$
            % \frac{\partial}{\partial t} A_0^{klj} = \underbrace{M_{A_0}^{klj} \left(\bar{A}_0^{klj}  - A_0^{klj} \right)/ \tau_0}_{F_\textrm{force}} + F_0^{klj} + F_\textrm{damp}^{klj}
            % $$
            %
            % - Topic: Computation
            % - Declaration: varargout = compute(wvt,varargin)
            % - Parameter A0bar: A0 'mean' value to relax to
            % - Parameter MA0: (optional) forcing mask, A0. 1s at the forced modes, 0s at the unforced modes. If it is left blank, then it will be produced using the nonzero values of A0bar
            % - Parameter tau0: (optional) relaxation time
            % - Returns varargout: cell array of returned variables
            arguments
                self WVNonlinearFluxQGForced {mustBeNonempty}
                A0bar (:,:,:) double {mustBeNonempty}
                options.MA0 (:,:,:) logical = abs(A0bar) > 1
                options.tau0 (1,1) double = 0
            end

            self.A0bar = A0bar;
            self.MA0 = options.MA0;
            self.tau0 = options.tau0;
        end

        function model_spectrum = setNarrowBandForcing(self, options)
            arguments
                self WVNonlinearFluxQGForced {mustBeNonempty}
                options.r (1,1) double
                options.k_r (1,1) double =(self.wvt.k(2)-self.wvt.k(1))*2
                options.k_f (1,1) double =(self.wvt.k(2)-self.wvt.k(1))*20
                options.j_f (1,1) double = 1
                options.u_rms (1,1) double = 0.2
                options.initialPV {mustBeMember(options.initialPV,{'none','narrow-band','full-spectrum'})} = 'narrow-band'
            end
            magicNumber = 0.0225;
            if isfield(options,"r")
                r = options.r;
                k_r = r/(magicNumber*options.u_rms);
            else
                r = magicNumber*options.u_rms*options.k_r; % 1/s bracket [0.02 0.025]
                k_r = options.k_r;
            end
            k_f = options.k_f;
            j_f = options.j_f;
            wvt = self.wvt;

            smallDampIndex = find(abs(self.damp(:,1)) > 1.1*abs(r),1,'first');
            fprintf('Small scale damping begins around k=%d dk. You have k_f=%d dk.\n',smallDampIndex-1,round(k_f/(wvt.k(2)-wvt.k(1))));

            
            deltaK = wvt.kRadial(2)-wvt.kRadial(1);
            self.MA0 = zeros(wvt.Nk,wvt.Nl,wvt.Nj);
            self.MA0(wvt.Kh > k_f-deltaK/2 & wvt.Kh < k_f+deltaK/2 & wvt.J == j_f) = 1;

            if strcmp(options.initialPV,'narrow-band') || strcmp(options.initialPV,'full-spectrum')
                u_rms = options.u_rms;

                m = 3/2; % We don't really know what this number is.
                kappa_epsilon = 0.5 * u_rms^2 / ( ((3*m+5)/(2*m+2))*k_r^(-2/3) - k_f^(-2/3) );
                model_viscous = @(k) kappa_epsilon * k_r^(-5/3 - m) * k.^m;
                model_inverse = @(k) kappa_epsilon * k.^(-5/3);
                model_forward = @(k) kappa_epsilon * k_f^(4/3) * k.^(-3);
                model_spectrum = @(k) model_viscous(k) .* (k<k_r) + model_inverse(k) .* (k >= k_r & k<=k_f) + model_forward(k) .* (k>k_f);

                % In this loop we set the energy level directly to wvt.A0
                % as computed by integrated the spectrum--so it is *not*
                % depth integrated and must be converted at the end.
                % ARand contains the random orientations
                kAxis = wvt.kRadial;
                dk = kAxis(2)-kAxis(1);
                ARand = wvt.generateHermitianRandomMatrix();
                for iK=1:(length(wvt.kRadial)-1)
                    indicesForK = find( kAxis(iK)-dk/2 <= wvt.Kh & wvt.Kh < kAxis(iK)+dk/2 & wvt.J == j_f  );
                    energy = integral(model_spectrum,max(kAxis(iK)-dk/2,0),kAxis(iK)+dk/2);
                    wvt.A0(indicesForK) = energy/length(indicesForK);
                    ARand(indicesForK) = ARand(indicesForK) /sqrt( sum(ARand(indicesForK) .* conj( ARand(indicesForK)))/length(indicesForK) );
                end
                wvt.A0 = WVTransform.makeHermitian(wvt.A0);

                AA = ~(wvt.maskForAliasedModes(jFraction=1));
                wvt.A0 = AA .* (sqrt(wvt.h .* wvt.A0) ./sqrt(wvt.A0_TE_factor)) .* ARand;
                wvt.A0(isnan(wvt.A0)) = 0;
                WVTransform.checkHermitian(wvt.A0);

                if strcmp(options.initialPV,'narrow-band')
                    wvt.A0 = (self.MA0) .* wvt.A0;
                else
                    fprintf('desired energy: %g, actual energy %g\n',0.5 * u_rms^2,wvt.geostrophicEnergy/wvt.h(j_f+1));
                end
            end
            self.A0bar = (self.MA0) .* wvt.A0;
            self.tau0 = 0;
        end

        function varargout = compute(self,wvt,varargin)
            varargout = cell(1,self.nVarOut);
            [varargout{:}] = compute@WVNonlinearFluxQG(self,wvt,varargin{:});

            if self.tau0 > 0
                varargout{1} = self.MA0.*(self.A0bar - wvt.A0)/self.tau0 + varargout{1};
            else
                varargout{1} = (~self.MA0) .* varargout{1};
            end
        end

        function writeToFile(self,ncfile,wvt)
            arguments
                self WVNonlinearFluxOperation {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            writeToFile@WVNonlinearFluxQG(self,ncfile,wvt);
            ncfile.addVariable('MA0',int8(self.MA0),{'k','l','j'});
            ncfile.addVariable('A0bar',self.A0bar,{'k','l','j'});
            ncfile.addVariable('tau0',wvt.tau0,{});
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
            nlFlux = WVNonlinearFluxQGForced(wvt,r=ncfile.attributes('r'),nu=ncfile.attributes('nu'),shouldUseBeta=(ncfile.attributes('beta')>0) );
            nlFlux.MA0 = logical(ncfile.readVariables('MA0'));
            nlFlux.A0bar = ncfile.readVariables('A0bar');
            nlFlux.tau0 = ncfile.readVariables('tau0');
        end
    end

end