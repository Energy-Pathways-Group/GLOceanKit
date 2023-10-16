classdef WVFluxQGUnforced < WVNonlinearFluxOperation
    properties
        PVA0 % conversion from A0 to PV
        A0PV % conversion from PV to A0
        beta
    end

    methods
        function self = WVFluxQGUnforced(wvt,options)
            arguments
                wvt WVTransform {mustBeNonempty}
                options.shouldUseBeta double {mustBeMember(options.shouldUseBeta,[0 1])} = 0 
                options.fluxName char = 'WVFluxQGUnforced'
                options.stateVariables WVVariableAnnotation = WVVariableAnnotation.empty()
            end
            fluxVar(1) = WVVariableAnnotation('F0',{'k','l','j'},'m/s', 'non-linear flux into A0');
            fluxVar(2) = WVVariableAnnotation('u',{'x','y','z'},'m/s', 'geostrophic velocity x-direction');
            fluxVar(2).attributes('standard_name') = 'eastward_sea_water_velocity';
            fluxVar(3) = WVVariableAnnotation('v',{'x','y','z'},'m/s', 'geostrophic velocity y-direction');
            fluxVar(3).attributes('standard_name') = 'northward_sea_water_velocity';
            fluxVar = cat(2,fluxVar,options.stateVariables);

            self@WVNonlinearFluxOperation(options.fluxName,fluxVar);
            self.doesFluxAp = 0;
            self.doesFluxAm = 0;
            self.doesFluxA0 = 1;
            
            AA = ~(wvt.maskForAliasedModes(jFraction=1));
            self.PVA0 = wvt.A0_QGPV_factor;
            self.A0PV = AA./self.PVA0;
            
            if options.shouldUseBeta == 1
                self.beta = 2 * 7.2921E-5 * cos( wvt.latitude*pi/180. ) / 6.371e6;
            else
                self.beta = 0;
            end
        end

        function varargout = compute(self,wvt,varargin)
            % Apply operator S---defined in (C4) in the manuscript
            Ubar = wvt.UA0 .* wvt.A0;
            Vbar = wvt.VA0 .* wvt.A0;
            PVbar = self.PVA0 .* wvt.A0;

            u_g = wvt.transformToSpatialDomainWithF(Ubar);
            v_g = wvt.transformToSpatialDomainWithF(Vbar);
%             QGPV = wvt.transformToSpatialDomainWithF(PVbar);
            PVx = wvt.transformToSpatialDomainWithF(sqrt(-1)*wvt.k.*PVbar);
            PVy = wvt.transformToSpatialDomainWithF(sqrt(-1)*shiftdim(wvt.l,-1).*PVbar);

            PVnl = u_g.*PVx + v_g.*(PVy+self.beta);
            F0 = -self.A0PV .* wvt.transformFromSpatialDomainWithF(PVnl);
            varargout = {F0,u_g,v_g};
        end

        function writeToFile(self,ncfile,wvt)
            arguments
                self WVNonlinearFluxOperation {mustBeNonempty}
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            ncfile.addAttribute('beta',self.beta)
        end

        function nlFlux = nonlinearFluxWithDoubleResolution(self,wvtX2)
            nlFlux = WVFluxQGUnforced(wvtX2,r=self.r,shouldUseBeta=(self.beta>0));
        end

        function flag = isequal(self,other)
            arguments
                self WVNonlinearFluxOperation
                other WVNonlinearFluxOperation
            end
            flag = isequal@WVNonlinearFluxOperation(self,other);
            flag = flag & isequal(self.PVA0, other.PVA0);
            flag = flag & isequal(self.beta, other.beta);
        end

    end

    methods (Static)
        function nlFlux = nonlinearFluxFromFile(ncfile,wvt)
            arguments
                ncfile NetCDFFile {mustBeNonempty}
                wvt WVTransform {mustBeNonempty}
            end
            nlFlux = WVFluxQGUnforced(wvt,shouldUseBeta=(ncfile.attributes('beta')>0) );
        end
    end

end