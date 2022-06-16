classdef SingleModeForcedDissipativeQGPVEMasked < SingleModeQGPVE
    properties
        k_r
        k_f
        u_rms

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
            end
            k_r = options.k_r;
            k_f = options.k_f;
            u_rms = options.u_rms;

            r = 0.1*(1/sqrt(5))*u_rms*k_r; % 1/s
            nu = (3/2)*(wvt.x(2)-wvt.x(1))*u_rms; % m^2/s

            self@SingleModeQGPVE(wvt,fluxName='SingleModeForcedDissipativeQGPVE',r=r,nu=nu,shouldUseBeta=options.shouldUseBeta);

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
            
            kappa_epsilon = 1/100;
            epsilon = u_rms^3 * (2*pi)^3 * k_r;
            desiredEnergyOld = kappa_epsilon * epsilon^(2/3) * k_f^(-5/3) * deltaK;

            kappa_epsilon = u_rms^2 / (5*k_r^(-2/3) - 2*k_f^(-2/3));
            desiredEnergy = kappa_epsilon * k_f^(-5/3) * deltaK;
            A0 = sqrt(wvt.h*desiredEnergy) * WaveVortexTransform.generateHermitianRandomMatrix(size(wvt.A0)) .* (F./sqrt(wvt.A0_TE_factor));

            wvt.A0 = A0;

            fprintf('desired energy: %g, actual energy %g. energy density: %g\n',desiredEnergy,wvt.geostrophicEnergy/wvt.h,desiredEnergy/deltaK);


            self.k_f = options.k_f;
            self.k_r = options.k_r;
            self.u_rms = options.u_rms;
        end

        function varargout = Compute(self,wvt,varargin)
            varargout = cell(1,self.nVarOut);
            [varargout{:}] = Compute@SingleModeQGPVE(self,wvt,varargin{:});
            F0 = varargout{1};
            F0 = self.EMA0 .* F0;
            varargout{1} = F0;
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

%         function initFromFile(self,ncfile,wvt)
%             arguments
%                 self NonlinearFluxOperation {mustBeNonempty}
%                 ncfile NetCDFFile {mustBeNonempty}
%                 wvt WaveVortexTransform {mustBeNonempty}
%             end
%             if isKey(ncfile.variableWithName,'beta')
%                 self.beta = ncfile.
%             end
%         end

    end

end