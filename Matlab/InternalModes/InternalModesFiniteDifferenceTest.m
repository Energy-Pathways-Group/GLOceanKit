classdef InternalModesFiniteDifferenceTest < InternalModesFiniteDifference
    % This class uses finite differencing of arbitrary order to compute the
    % internal wave modes. See InternalModesBase for basic usage
    % information.
    %
    %   The class takes the name/value pair 'orderOfAccuracy' (default
    %   value of 4) in order to set the order of accuracy of the finite
    %   differencing matrix. The matrix is constructed using the weights
    %   algorithm described by Bengt Fornberg in 'Calaculation of weight in
    %   finite difference formulas', SIAM review, 1998.
    %
    %   Setting the orderOfAccuracy does tended to improve the quality of
    %   the solution, but does tend to have strange effects when the order
    %   gets high relative to the number of grid points.
    %
    %   If you request a different output grid than input grid, the
    %   solutions are mapped to the output grid with spline interpolation.
    %
    %   See also INTERNALMODES, INTERNALMODESSPECTRAL,
    %   INTERNALMODESDENSITYSPECTRAL, INTERNALMODESWKBSPECTRAL, and
    %   INTERNALMODESBASE.
    %
    %
    %   Jeffrey J. Early
    %   jeffrey@jeffreyearly.com
    %
    %   March 14th, 2017        Version 1.0
    
    methods
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        function self = InternalModesFiniteDifferenceTest(rho, z_in, z_out, latitude, varargin)
            % Initialize with either a grid or analytical profile.
            self@InternalModesFiniteDifference(rho,z_in,z_out,latitude,varargin{:});
        end
        
        function [F,G,h] = ModesAtFrequency(self, omega )
            % Return the normal modes and eigenvalue at a given frequency.
            
            A = self.Diff2;
            B = -diag(self.N2_z_diff - omega*omega)/self.g;
            
            a = self.N2_z_diff - omega*omega;
            a(a>=0) = 1;
            a(a<0) = 0;
            turningIndex = find(diff(a)~=0);
            A(turningIndex,:) = 0;
            B(turningIndex,:) = 0;
            A(turningIndex,turningIndex-1) = 1;
            A(turningIndex,turningIndex) = -1;
            
            % Bottom boundary condition (always taken to be G=0)
            % NOTE: we already chose the correct BCs when creating the
            % Diff2 matrix
            A(1,:) = self.Diff2(1,:);
            B(1,:) = 0;
            
            % Surface boundary condition
            if strcmp(self.upperBoundary, 'free_surface')
                % G_z = \frac{1}{h_j} G at the surface
                B(end,end)=1;
            elseif strcmp(self.upperBoundary, 'rigid_lid')
                % G=0 at the surface (note we chose this BC when creating Diff2)
                A(end,:) = self.Diff2(end,:);
                B(end,end)=0;
            end
            
            h_func = @(lambda) 1.0 ./ lambda;
            [F,G,h] = ModesFromGEP(self,A,B,h_func);
        end
        
    end
end

