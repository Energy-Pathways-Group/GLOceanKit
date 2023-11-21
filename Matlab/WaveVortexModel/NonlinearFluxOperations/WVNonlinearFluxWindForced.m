classdef WVNonlinearFluxWindForced < WVNonlinearFluxForced

    properties
        f_x, f_y
        slabDampTime

        u_wind
        v_wind
        t_wind
        maxT
        slabDepth
    end

    methods
        function self = WVNonlinearFluxWindForced(wvt,options,newOptions)
            arguments
                wvt WVTransform {mustBeNonempty}
                options.uv_damp (1,1) double 
                options.w_damp (1,1) double % characteristic speed used to set the damping. Try using wMax
                options.nu_xy (1,1) double
                options.nu_z (1,1) double
                options.r (1,1) double {mustBeNonnegative} = 0 % linear bottom friction, try 1/(200*86400) https://www.nemo-ocean.eu/doc/node70.html
                options.shouldAntialias double = 1

                newOptions.slabDampTime double = 10*86400
            end
            forcedArgs = namedargs2cell(options);
            self@WVNonlinearFluxForced(wvt,forcedArgs{:});
            
            % Allow energy fluxes at all modes
            self.slabDampTime = newOptions.slabDampTime;

            wind = load('wind');
            nan_indices = isnan(wind.u) | isnan(wind.v);
            self.u_wind = wind.u;
            self.v_wind = wind.v;
            self.u_wind(nan_indices) = 0;
            self.v_wind(nan_indices) = 0;
            self.t_wind = wind.t;
            self.maxT = max(wind.t);
            self.slabDepth = abs(wvt.z(end-1));

            fprintf('WVNonlinearFluxWindForced: effective slab depth: %.1f m, decay time: %.1f days\n',self.slabDepth,self.slabDampTime/86400);
        end

        function varargout = spatialFlux(self,wvt,varargin)
            [uNL,vNL,nNL,phase,uvw] = spatialFlux@WVNonlinearFluxForced(self,wvt,varargin{:});
            
            indices = wvt.Nz;
            uSurface = uvw.u(:,:,indices);
            vSurface = uvw.v(:,:,indices);
            [u_force,v_force] = self.slabForcingFromWinds(wvt.t,uSurface,vSurface);
            uNL(:,:,indices) = uNL(:,:,indices) + u_force - uSurface / self.slabDampTime;
            vNL(:,:,indices) = vNL(:,:,indices) + v_force - vSurface / self.slabDampTime;

            if nargout == 5
                varargout = {uNL,vNL,nNL,phase,uvw};
            else
                varargout = {uNL,vNL,nNL,phase};
            end
        end

        function [u_force,v_force] = slabForcingFromWinds(self,t,u_ocean,v_ocean)
            % The wind stress is computed from the wind vector using the empirical coefficients
            % found on page 1520 of Yelland, et al., 1998, JPO. Technically this is only valid
            % for speeds between 6 and 25 m/s.
            u = interp1(self.t_wind,self.u_wind,mod(t,self.maxT)) - u_ocean;
            v = interp1(self.t_wind,self.v_wind,mod(t,self.maxT)) - v_ocean;
            rho_air = 1.25;
            speed_wind = sqrt(u.^2 + v.^2);
                
            % equation 3
            drag_coefficient = (1e-3)*(0.50 + 0.071*speed_wind);

            % stress tau_x, tau_y returned in units of N/m^2, or kg m^{-1} s^{-2}
            % divide by rho0*slabDepth convert stress to forcing, m s^{-2}
            u_force = rho_air * drag_coefficient .* speed_wind .* u/(self.wvt.rho0 * self.slabDepth); 
            v_force = rho_air * drag_coefficient .* speed_wind .* v/(self.wvt.rho0 * self.slabDepth); 
        end
    end
end