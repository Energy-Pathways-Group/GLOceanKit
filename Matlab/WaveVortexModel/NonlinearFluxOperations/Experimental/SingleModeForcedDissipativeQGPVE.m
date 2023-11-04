classdef SingleModeForcedDissipativeQGPVE < SingleModeQGPVE
    properties
        k_f
        f_zeta
        F
    end

    methods
        function self = SingleModeForcedDissipativeQGPVE(wvt,options)
            arguments
                wvt WVTransform {mustBeNonempty}
                options.shouldUseBeta double {mustBeMember(options.shouldUseBeta,[0 1])} = 0 
                options.k_f (1,1) double
                options.f_zeta (1,1) double
                options.r (1,1) double
                options.nu (1,1) double
            end

            [K,L,~] = ndgrid(wvt.k,wvt.l,wvt.j);
            Kh = sqrt(K.*K + L.*L);
            deltaK = wvt.kRadial(2)-wvt.kRadial(1);

            % This creates an annulus with unit variance.
            Fzeta = zeros(size(Kh));
            Fzeta(Kh > options.k_f-deltaK & Kh < options.k_f+deltaK) = 1;
            Fzeta = WVTransform.makeHermitian(Fzeta);
            Fzeta = Fzeta/sqrt(sum(sum(Fzeta.*conj(Fzeta))));

            self@SingleModeQGPVE(wvt,fluxName='SingleModeForcedDissipativeQGPVE',r=options.r,nu=options.nu,shouldUseBeta=options.shouldUseBeta);

            self.k_f = options.k_f;
            self.f_zeta = options.f_zeta;
            self.F = self.f_zeta*(wvt.f*wvt.h ./(wvt.Omega .* wvt.Omega)).*Fzeta;
        end

        function varargout = compute(self,wvt,varargin)
            varargout = cell(1,self.nVarOut);
            [varargout{:}] = compute@SingleModeQGPVE(self,wvt,varargin{:});
            F0 = varargout{1};
            F0 = F0 + wvt.generateHermitianRandomMatrix().*self.F;
            varargout{1} = F0;
        end

        function writeToFile(self,ncfile,wvt)
            arguments
                self WVNonlinearFluxOperation {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end

            writeToFile@SingleModeQGPVE(self,ncfile,wvt);
            ncfile.addAttribute('k_f',self.k_f)
            ncfile.addAttribute('f_zeta',self.f_zeta)
        end

%         function initFromFile(self,ncfile,wvt)
%             arguments
%                 self WVNonlinearFluxOperation {mustBeNonempty}
%                 ncfile NetCDFFile {mustBeNonempty}
%                 wvt WVTransform {mustBeNonempty}
%             end
%             if isKey(ncfile.variableWithName,'beta')
%                 self.beta = ncfile.
%             end
%         end

    end

end