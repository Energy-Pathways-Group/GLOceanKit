classdef SingleModeForcedDissipativeQGPVEMasked < SingleModeQGPVE
    properties
        k_f
        EMA0

        model_spectrum
        k_r
    end

    methods
        function self = SingleModeForcedDissipativeQGPVEMasked(wvt,options)
            arguments
                wvt WVTransform {mustBeNonempty}
                options.r (1,1) double
                options.nu (1,1) double
                options.k_r (1,1) double =(wvt.k(2)-wvt.k(1))*2
                options.k_f (1,1) double =(wvt.k(2)-wvt.k(1))*20
                options.u_rms (1,1) double = 0.2
                options.shouldUseBeta double {mustBeMember(options.shouldUseBeta,[0 1])} = 0
                options.initialPV {mustBeMember(options.initialPV,{'none','narrow-band','full-spectrum'})} = 'narrow-band'
            end

            if isfield(options,"nu")
                nu = options.nu;
            else
                % The correct setting for nu would be the following,
                % nu = (3/2)*(wvt.x(2)-wvt.x(1))*options.uMax/(pi^2);
                % but we don't know what u_max is, do we?
                nu = (3/2)*(wvt.x(2)-wvt.x(1))*options.u_rms;
            end

            magicNumber = 0.0225;
            if isfield(options,"r")
                r = options.r;
                k_r = r/(magicNumber*options.u_rms);
            else
                r = magicNumber*options.u_rms*options.k_r; % 1/s bracket [0.02 0.025]
                k_r = options.k_r;
            end

            fluxVar(1) = WVVariableAnnotation('F0_psi',{'k','l','j'},'m/s', 'forcing function applied to the vortex coefficients',isComplex=1);
            fluxVar(2) = WVVariableAnnotation('F_psi',{'x','y','z'},'1/s^2', 'forcing function applied to the QGPVE written in terms of psi');
            self@SingleModeQGPVE(wvt,fluxName='SingleModeForcedDissipativeQGPVE',r=r,nu=nu,shouldUseBeta=options.shouldUseBeta,stateVariables=fluxVar);
            self.k_f = options.k_f;
            self.k_r = k_r;

            smallDampIndex = find(abs(self.damp(:,1)) > 1.1*abs(r),1,'first');
            fprintf('Small scale damping begins around k=%d dk. You have k_f=%d dk.\n',smallDampIndex-1,round(self.k_f/(wvt.k(2)-wvt.k(1))));

            % Create an energy mask, so we don't remove energy from the
            % annulus.
            deltaK = wvt.kRadial(2)-wvt.kRadial(1);
            self.EMA0 = ones(wvt.Nk,wvt.Nl,wvt.Nj);
            self.EMA0(wvt.Kh > self.k_f-deltaK/2 & wvt.Kh < self.k_f+deltaK/2) = 0;
            
            if strcmp(options.initialPV,'narrow-band') || strcmp(options.initialPV,'full-spectrum')
                k_f = self.k_f;
                u_rms = options.u_rms;
                
                m = 3/2; % We don't really know what this number is.
                kappa_epsilon = 0.5 * u_rms^2 / ( ((3*m+5)/(2*m+2))*k_r^(-2/3) - k_f^(-2/3) );
                model_viscous = @(k) kappa_epsilon * k_r^(-5/3 - m) * k.^m;
                model_inverse = @(k) kappa_epsilon * k.^(-5/3);
                model_forward = @(k) kappa_epsilon * k_f^(4/3) * k.^(-3);
                self.model_spectrum = @(k) model_viscous(k) .* (k<k_r) + model_inverse(k) .* (k >= k_r & k<=k_f) + model_forward(k) .* (k>k_f);

                kAxis = wvt.kRadial;
                dk = kAxis(2)-kAxis(1);
                ARand = wvt.generateHermitianRandomMatrix();
                for iK=1:(length(wvt.kRadial)-1)
                    indicesForK = find( kAxis(iK)-dk/2 <= wvt.Kh & wvt.Kh < kAxis(iK)+dk/2   );
                    energy = integral(self.model_spectrum,max(kAxis(iK)-dk/2,0),kAxis(iK)+dk/2);
                    wvt.A0(indicesForK) = energy/length(indicesForK);
                    ARand(indicesForK) = ARand(indicesForK) /sqrt( sum(ARand(indicesForK) .* conj( ARand(indicesForK)))/length(indicesForK) );
                end
                wvt.A0 = WVTransform.makeHermitian(wvt.A0);

                AA = ~(wvt.maskForAliasedModes(jFraction=1));
                wvt.A0 = AA .* (sqrt(wvt.h * wvt.A0) ./sqrt(wvt.A0_TE_factor)) .* ARand;
                WVTransform.checkHermitian(wvt.A0);

                if strcmp(options.initialPV,'narrow-band')
                    wvt.A0 = (~self.EMA0) .* wvt.A0;
                else
                    fprintf('desired energy: %g, actual energy %g\n',0.5 * u_rms^2,wvt.geostrophicEnergy/wvt.h);
                end
            end

            
        end
        

        function varargout = compute(self,wvt,varargin)
            varargout = cell(1,self.nVarOut-2);
            [varargout{:}] = compute@SingleModeQGPVE(self,wvt,varargin{:});
            F0 = varargout{1};
            F0_psi = (~self.EMA0) .* F0;
            forcing = wvt.transformToSpatialDomainWithF(self.PVA0 .* F0_psi);
            F0 = self.EMA0 .* F0;
            varargout{1} = F0;
            varargout{end+1} = F0_psi;
            varargout{end+1} = forcing;
        end

        function writeToFile(self,ncfile,wvt)
            arguments
                self WVNonlinearFluxOperation {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end

            writeToFile@SingleModeQGPVE(self,ncfile,wvt);
            ncfile.addAttribute('k_f',self.k_f)
            ncfile.addAttribute('k_r',self.k_r)
        end

        function nlFlux = nonlinearFluxWithResolutionOfTransform(self,wvtX2)
            nlFlux = SingleModeForcedDissipativeQGPVEMasked(wvtX2,initialPV='none', k_f=self.k_f,r=self.r,nu=self.nu/2);
        end

        function flag = isequal(self,other)
            arguments
                self WVNonlinearFluxOperation
                other WVNonlinearFluxOperation
            end
            flag = isequal@SingleModeQGPVE(self,other);
            flag = flag & isequal(self.k_f,other.k_f);
            flag = flag & isequal(self.EMA0,other.EMA0);
        end
    end

    methods (Static)
        function nlFlux = nonlinearFluxFromFile(ncfile,wvt)
            arguments
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            nlFlux = SingleModeForcedDissipativeQGPVEMasked(wvt,initialPV='none', k_f=ncfile.attributes('k_f'),r=ncfile.attributes('r'),nu=ncfile.attributes('nu') );
        end
    end

end