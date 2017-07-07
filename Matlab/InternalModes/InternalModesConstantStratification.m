classdef InternalModesConstantStratification < InternalModesBase
    properties (Access = public)
        N0
    end
    
    methods
        function self = InternalModesConstantStratification(rho, z_in, z_out, latitude, varargin) 
            self@InternalModesBase(rho,z_in,z_out,latitude, varargin{:});
            self.N0 = InternalModesConstantStratification.BuoyancyFrequencyFromConstantStratification(rho,z_in,self.rho0,self.g);
        end
    end
    
    methods (Static)
        function flag = IsStratificationConstant(rho,z_in)
            if isa(rho,'function_handle') == true
                if numel(z_in) ~= 2
                    error('When using a function handle, z_domain must be an array with two values: z_domain = [z_bottom z_surface];')
                end
                rho0 = rho(max(z_in));
                max_ddrho = max(abs(diff(diff(rho(linspace(min(z_in),max(z_in),100))/rho0 ))));
            elseif isa(rho,'numeric') == true
                rho0 = max(rho);
                max_ddrho = max(abs(diff(diff( rho/rho0 ))));
            end
            
            flag = max_ddrho < 1e-7;
        end
        
        function N0 = BuoyancyFrequencyFromConstantStratification(rho,z_in,rho0,g)
            if isa(rho,'function_handle') == true
                drhodz = (rho(max(z_in)) - rho(min(z_in)))/( max(z_in) - min(z_in) );
            elseif isa(rho,'numeric') == true
                drhodz = (min(rho)-max(rho))/(max(z_in)-min(z_in));
            end
            N0 = sqrt(-g*drhodz/rho0);
        end
        
    end
end