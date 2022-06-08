classdef SingleModeQGPVE < NonlinearFluxOperation
    properties
        PVA0 % conversion from A0 to PV
        A0PV % conversion from PV to A0
        beta
    end

    methods
        function self = SingleModeQGPVE(wvt,options)
            arguments
                wvt WaveVortexTransform {mustBeNonempty}
                options.shouldUseBeta double {mustBeMember(options.shouldUseBeta,[0 1])} = 0 
            end
            fluxVar(1) = StateVariable('F0',{'k','l','j'},'m/s', 'non-linear flux into A0');
            fluxVar(2) = StateVariable('u_g',{'x','y','z'},'m/s', 'geostrophic velocity x-direction');
            fluxVar(3) = StateVariable('v_g',{'x','y','z'},'m/s', 'geostrophic velocity y-direction');

            self@NonlinearFluxOperation('SingleModeQGPVE',fluxVar);
            self.doesFluxAp = 0;
            self.doesFluxAm = 0;
            self.doesFluxA0 = 1;

            self.PVA0 = - wvt.PP .* wvt.Omega .* wvt.Omega / (wvt.h * wvt.f0);
            self.A0PV = 1./self.PVA0;
            if options.shouldUseBeta == 1
                self.beta = 2 * 7.2921E-5 * cos( wvt.latitude*pi/180. ) / 6.371e6;
            else
                self.beta = 0;
            end
        end

        function varargout = Compute(self,wvt,varargin)
            % Apply operator S---defined in (C4) in the manuscript
            Ubar = wvt.UA0 .* wvt.A0;
            Vbar = wvt.VA0 .* wvt.A0;
            PVbar = self.PVA0 .* wvt.A0;

            u_g = ifft(ifft(Ubar,wvt.Nx,1),wvt.Ny,2,'symmetric');
            v_g = ifft(ifft(Vbar,wvt.Nx,1),wvt.Ny,2,'symmetric');
            PVx = ifft( sqrt(-1)*wvt.k.*ifft(PVbar,wvt.Ny,2), wvt.Nx, 1,'symmetric');
            PVy = ifft( sqrt(-1)*shiftdim(wvt.l,-1).*ifft(PVbar,wvt.Nx,1), wvt.Ny, 2,'symmetric');

            PVnl = u_g.*PVx + v_g.*(PVy+self.beta);
            F0 = -self.A0PV .* wvt.TransformFromSpatialDomainWithF(PVnl);
            varargout = {F0,u_g,v_g};
        end

        function writeToFile(self,ncfile,wvt)
            arguments
                self NonlinearFluxOperation {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                wvt WaveVortexTransform {mustBeNonempty}
            end
            ncfile.addAttribute('beta',self.beta)
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