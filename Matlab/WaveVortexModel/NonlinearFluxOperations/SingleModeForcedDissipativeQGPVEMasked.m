classdef SingleModeForcedDissipativeQGPVEMasked < SingleModeQGPVE
    properties
        k_r
        k_f
        u_rms

        model_spectrum

        EMA0
    end

    methods
        function self = SingleModeForcedDissipativeQGPVEMasked(wvt,options)
            arguments
                wvt WaveVortexTransform {mustBeNonempty}
                options.k_r (1,1) double =(wvt.k(2)-wvt.k(1))*2
                options.k_f (1,1) double =(wvt.k(2)-wvt.k(1))*20
                options.u_rms (1,1) double = 0.2
                options.shouldUseBeta double {mustBeMember(options.shouldUseBeta,[0 1])} = 0
                options.initialPV {mustBeMember(options.initialPV,{'none','narrow-band','full-spectrum'})} = 'narrow-band'
            end
            k_r = options.k_r;
            k_f = options.k_f;
            u_rms = options.u_rms;
            
            % experimenting with the prefactor c_r in r=c_r*(1/sqrt(5))*u_rms*k_r
            % I find 0.1 *over* damped, but 0.05 underdamped, at least by a
            % little.
            % Switching to 256 changes this... needed to free-up more
            % inertial range. 0.05 runs great
            r = 0.04*(1/sqrt(5))*u_rms*k_r; % 1/s
            r = 0.0225*u_rms*k_r; % 1/s bracket [0.02 0.025]
%             r = 0.04*u_rms*k_r; % 1/s bracket [0.02 0.025]
            nu = (3/2)*(wvt.x(2)-wvt.x(1))*u_rms; % m^2/s

            fluxVar(1) = StateVariable('F0_psi',{'k','l','j'},'m/s', 'forcing function applied to the vortex coefficients',isComplex=1);
            fluxVar(2) = StateVariable('F_psi',{'x','y','z'},'1/s^2', 'forcing function applied to the QGPVE written in terms of psi');

            self@SingleModeQGPVE(wvt,fluxName='SingleModeForcedDissipativeQGPVE',r=r,nu=nu,shouldUseBeta=options.shouldUseBeta,stateVariables=fluxVar);
            smallDampIndex = find(abs(self.damp(:,1)) > 1.1*abs(r),1,'first');
            fprintf('Small scale damping begins around k=%d dk. You have k_f=%d dk.\n',smallDampIndex-1,round(k_f/(wvt.k(2)-wvt.k(1))));

            [K,L,~] = ndgrid(wvt.k,wvt.l,wvt.j);
            Kh = sqrt(K.*K + L.*L);
            deltaK = wvt.kRadial(2)-wvt.kRadial(1);

            % This creates an annulus with unit variance.
            F = zeros(size(Kh));
            F(Kh > options.k_f-deltaK/2 & Kh < options.k_f+deltaK/2) = 1;
            F = WaveVortexTransform.makeHermitian(F);
            F = F/sqrt(sum(sum(F.*conj(F))));

            % Create an energy mask, so we don't remove energy from the
            % annulus.
            self.EMA0 = ones(wvt.Nk,wvt.Nl,wvt.Nj);
            self.EMA0(F>0) = 0;
            
            if strcmp(options.initialPV,'narrow-band') || strcmp(options.initialPV,'full-spectrum')
                m = 3/2; % We don't really know what this number is.
                kappa_epsilon = 0.5 * u_rms^2 / ( ((3*m+5)/(2*m+2))*k_r^(-2/3) - k_f^(-2/3) );
                model_viscous = @(k) kappa_epsilon * k_r^(-5/3 - m) * k.^m;
                model_inverse = @(k) kappa_epsilon * k.^(-5/3);
                model_forward = @(k) kappa_epsilon * k_f^(4/3) * k.^(-3);
                self.model_spectrum = @(k) model_viscous(k) .* (k<k_r) + model_inverse(k) .* (k >= k_r & k<=k_f) + model_forward(k) .* (k>k_f);

                kAxis = wvt.kRadial;
                dk = kAxis(2)-kAxis(1);
                ARand = WaveVortexTransform.generateHermitianRandomMatrix(size(wvt.A0));
                for iK=1:(length(wvt.kRadial)-1)
                    indicesForK = find( kAxis(iK)-dk/2 <= squeeze(Kh(:,:,1)) & squeeze(Kh(:,:,1)) < kAxis(iK)+dk/2   );
                    energy = integral(self.model_spectrum,max(kAxis(iK)-dk/2,0),kAxis(iK)+dk/2);
                    wvt.A0(indicesForK) = energy/length(indicesForK);
                    ARand(indicesForK) = ARand(indicesForK) /sqrt( sum(ARand(indicesForK) .* conj( ARand(indicesForK)))/length(indicesForK) );
                end
                wvt.A0 = WaveVortexTransform.makeHermitian(wvt.A0);

                %                     dk = wvt.k(2)-wvt.k(1);
                AA = ~(wvt.MaskForAliasedModes(jFraction=1));
                wvt.A0 = AA .* (sqrt(wvt.h * wvt.A0) ./sqrt(wvt.A0_TE_factor)) .* ARand;

                if strcmp(options.initialPV,'narrow-band')
                    wvt.A0(~(F>0)) = 0;
                else
                    fprintf('desired energy: %g, actual energy %g\n',0.5 * u_rms^2,wvt.geostrophicEnergy/wvt.h);
                end
        
%                 else
%                     kappa_epsilon = 1/100;
%                     epsilon = u_rms^3 * (2*pi)^3 * k_r;
%                     desiredEnergyOld = kappa_epsilon * epsilon^(2/3) * k_f^(-5/3) * deltaK;
% 
%                     kappa_epsilon = u_rms^2 / (5*k_r^(-2/3) - 2*k_f^(-2/3));
%                     desiredEnergy = kappa_epsilon * k_f^(-5/3) * deltaK;
%                     A0 = sqrt(wvt.h*desiredEnergy) * WaveVortexTransform.generateHermitianRandomMatrix(size(wvt.A0)) .* (F./sqrt(wvt.A0_TE_factor));
% 
%                     wvt.A0 = A0;
% 
%                     fprintf('desired energy: %g, actual energy %g. energy density: %g\n',desiredEnergy,wvt.geostrophicEnergy/wvt.h,desiredEnergy/deltaK);
%                 end
            end

            self.k_f = options.k_f;
            self.k_r = options.k_r;
            self.u_rms = options.u_rms;
        end

        function varargout = Compute(self,wvt,varargin)
            varargout = cell(1,self.nVarOut-2);
            [varargout{:}] = Compute@SingleModeQGPVE(self,wvt,varargin{:});
            F0 = varargout{1};
            F0_psi = (~self.EMA0) .* F0;
            forcing = ifft(ifft( self.PVA0 .* F0_psi,wvt.Nx,1),wvt.Ny,2,'symmetric');
            F0 = self.EMA0 .* F0;
            varargout{1} = F0;
            varargout{end+1} = F0_psi;
            varargout{end+1} = forcing;
        end

        function writeToFile(self,ncfile,wvt)
            arguments
                self NonlinearFluxOperation {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                wvt WaveVortexTransform {mustBeNonempty}
            end

            writeToFile@SingleModeQGPVE(self,ncfile,wvt);
            ncfile.addAttribute('k_f',self.k_f)
            ncfile.addAttribute('k_r',self.k_r)
            ncfile.addAttribute('u_rms',self.u_rms)
        end

        function nlFlux = nonlinearFluxWithDoubleResolution(self,wvtX2)
            nlFlux = SingleModeForcedDissipativeQGPVEMasked(wvtX2,initialPV='none', k_f=self.k_f,k_r=self.k_r,u_rms=self.u_rms);
        end
    end

    methods (Static)
        function nlFlux = nonlinearFluxFromFile(ncfile,wvt)
            arguments
                ncfile NetCDFFile {mustBeNonempty}
                wvt WaveVortexTransform {mustBeNonempty}
            end
            nlFlux = SingleModeForcedDissipativeQGPVEMasked(wvt,initialPV='none', k_f=ncfile.attributes('k_f'),k_r=ncfile.attributes('k_r'),u_rms=ncfile.attributes('u_rms') );
        end
    end

end