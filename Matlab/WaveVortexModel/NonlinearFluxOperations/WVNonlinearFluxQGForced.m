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
    % To initialize the WVNonlinearFluxQGForced,
    %
    % ```matlab
    % model = WVModel(wvt,nonlinearFlux=WVNonlinearFluxQGForced(wvt,shouldUseBeta=1,uv_damp=wvt.uMax));
    % ```
    %
    % - Topic: Initializing
    % - Declaration: WVNonlinearFluxQGForced < [WVNonlinearFluxQG](/classes/wvnonlinearfluxqg/)
    properties
        MA0         % Forcing mask, A0. 1s at the forced modes, 0s at the unforced modes
        A0bar = []  % A0 'mean' value to relax to
        tau0        % relaxation time
    end

    methods
        function self = WVNonlinearFluxQGForced(wvt,options,newOptions)
            % initialize 3D quasigeostrophic potential vorticity flux
            %
            % - Declaration: nlFlux = WVNonlinearFluxQGForced(wvt,options)
            % - Parameter wvt: a WVTransform instance
            % - Parameter shouldUseBeta: (optional) a Boolean indicating whether or not to include beta in the flux
            % - Parameter uv_damp: (optional) characteristic speed used to set the damping. Try using wvt.uMax
            % - Parameter r: (optional) bottom friction
            % - Parameter nu_xy: (optional) coefficient for damping
            % - Returns nlFlux: a QGPVE instance
            arguments
                wvt WVTransform {mustBeNonempty}
                options.shouldUseBeta double {mustBeMember(options.shouldUseBeta,[0 1])} = 0 
                options.uv_damp (1,1) double = 0.25 % characteristic speed used to set the damping. Try using uMax.
                options.r (1,1) double = 0
                options.fluxName char = 'WVNonlinearFluxQGForced'
                options.nu_xy (1,1) double
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
                options.MA0 (:,:,:) logical = abs(A0bar) > 0
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
                options.u_rms (1,1) double = 0.2 % set the *total* energy (not just kinetic) equal to 0.5*u_rms^2
                options.initialPV {mustBeMember(options.initialPV,{'none','narrow-band','full-spectrum'})} = 'narrow-band'
            end
            
            if ~self.wvt.isBarotropic
                % the idea is to set the energy at the sea-surface and
                % so we need to know the relative amplitude of this
                % mode at the surface.
                F = self.wvt.FinvMatrix;
                surfaceMag = 1/F(end,options.j_f+1);
                sbRatio = abs(F(end,options.j_f+1)/F(1,options.j_f+1));
                sbRatio = 1; % should we change the damping scale? Or no?
                h = self.wvt.h(options.j_f+1);
            else
                surfaceMag = 1;
                sbRatio = 1;
                h = self.wvt.h;
                options.j_f = 0;
            end

            magicNumber = 0.0225;
            if isfield(options,"r")
                self.r = options.r;
                k_r = self.r/(magicNumber*options.u_rms);
            else
                r = magicNumber*sbRatio*options.u_rms*options.k_r; % 1/s bracket [0.02 0.025]
                fprintf('1/r is %.1f days, switching to %.1f days\n',1/(self.r*86400),1/(r*86400));
                self.r = r;
                k_r = options.k_r;
            end
            k_f = options.k_f;
            j_f = options.j_f;
            wvt = self.wvt;

            smallDampIndex = find(abs(self.damp(:,1)) > 1.1*abs(self.r),1,'first');
            fprintf('(k_r=%.2f dk, k_f=%d dk, k_nu=%d dk.\n',k_r/(wvt.k(2)-wvt.k(1)),round(k_f/(wvt.k(2)-wvt.k(1))),smallDampIndex-1);
            % fprintf('Small scale damping begins around k=%d dk. You have k_f=%d dk.\n',smallDampIndex-1,round(k_f/(wvt.k(2)-wvt.k(1))));

            
            deltaK = wvt.kRadial(2)-wvt.kRadial(1);
            self.MA0 = zeros(wvt.Nk,wvt.Nl,wvt.Nj);
            self.MA0(wvt.Kh > k_f-deltaK/2 & wvt.Kh < k_f+deltaK/2 & wvt.J == j_f) = 1;

            if strcmp(options.initialPV,'narrow-band') || strcmp(options.initialPV,'full-spectrum')
                u_rms = surfaceMag * options.u_rms;

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
                    u = wvt.seaSurfaceU;
                    v = wvt.seaSurfaceV;
                    zeta = wvt.seaSurfaceHeight;
                    KE = mean(mean(0.5*(u.^2+v.^2)));
                    PE = mean(mean(0.5*(9.81*zeta.^2)/h));
                    u_rms_surface = mean(mean(sqrt(u.^2+v.^2)));
                    fprintf("surface u_rms: %.2g cm/s\n",100*u_rms_surface);
                    fprintf("surface energy, %g.\n",KE+PE);
                    fprintf('desired energy: %g, actual energy %g\n',0.5 * u_rms^2,wvt.geostrophicEnergy/h);
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
            ncfile.addVariable('tau0',self.tau0,{});
        end

        function nlFlux = nonlinearFluxWithResolutionOfTransform(self,wvtX2)
            ratio = wvtX2.Nk/self.wvt.Nk;
            nlFlux = WVNonlinearFluxQGForced(wvtX2,r=self.r,shouldUseBeta=(self.beta>0),nu_xy=self.nu_xy/ratio);
            if ~isempty(self.MA0)
                nlFlux.MA0 = WVTransform.spectralVariableWithResolution(self.MA0,[wvtX2.Nk wvtX2.Nl wvtX2.Nj]);
            end
            if ~isempty(self.A0bar)
                nlFlux.A0bar = WVTransform.spectralVariableWithResolution(self.A0bar,[wvtX2.Nk wvtX2.Nl wvtX2.Nj]);
            end
            nlFlux.tau0 = self.tau0;
        end

        function nlFlux = nonlinearFluxWithResolutionForTransform(self,wvtX2)
            ratio = wvtX2.Nk/self.wvt.Nk;
            nlFlux = WVNonlinearFluxQGForced(wvtX2,r=self.r,shouldUseBeta=(self.beta>0),nu_xy=self.nu_xy/ratio);
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
            nlFlux = WVNonlinearFluxQGForced(wvt,r=ncfile.attributes('r'),nu_xy=ncfile.attributes('nu_xy'),shouldUseBeta=(ncfile.attributes('beta')>0) );
            nlFlux.MA0 = logical(ncfile.readVariables('MA0'));
            nlFlux.A0bar = ncfile.readVariables('A0bar');
            nlFlux.tau0 = ncfile.readVariables('tau0');
        end
    end

end