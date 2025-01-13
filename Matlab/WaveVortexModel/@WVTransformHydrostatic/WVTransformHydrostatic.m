classdef WVTransformHydrostatic < WVTransformAbstractHydrostatic & WVWaveComponents & WVInertialOscillationMethods & WVGeostrophicMethods & WVMeanDensityAnomalyMethods & WVInternalGravityWaveMethods
    % A class for disentangling hydrostatic waves and vortices in variable stratification
    %
    % To initialization an instance of the WVTransformHydrostatic class you
    % must specific the domain size, the number of grid points and *either*
    % the density profile or the stratification profile.
    %
    % ```matlab
    % N0 = 3*2*pi/3600;
    % L_gm = 1300;
    % N2 = @(z) N0*N0*exp(2*z/L_gm);
    % wvt = WVTransformHydrostatic([100e3, 100e3, 4000],[64, 64, 65], N2=N2,latitude=30);
    % ```
    %
    % - Topic: Initialization
    % - Topic: Primary flow components
    % - Topic: Stratification
    % - Topic: Stratification — Vertical modes
    % - Topic: Stratification — Validation
    % - Topic: Initial conditions
    % - Topic: Energetics of flow components
    % - Topic: Operations
    %
    % - Declaration: classdef WVTransformHydrostatic < [WVTransform](/classes/wvtransform/)
    properties (Dependent)
        h_0  % [Nj 1]
        h_pm  % [Nj 1]
    end

    methods
        function self = WVTransformHydrostatic(Lxyz, Nxyz, options)
            % create a wave-vortex transform for variable stratification
            %
            % Creates a new instance of the WVTransformHydrostatic class
            % appropriate for disentangling hydrostatic waves and vortices
            % in variable stratification
            %
            % You must initialization by passing *either* the density
            % profile or the stratification profile.
            %
            % - Topic: Initialization
            % - Declaration: wvt = WVTransformHydrostatic(Lxyz, Nxyz, options)
            % - Parameter Lxyz: length of the domain (in meters) in the three coordinate directions, e.g. [Lx Ly Lz]
            % - Parameter Nxyz: number of grid points in the three coordinate directions, e.g. [Nx Ny Nz]
            % - Parameter rho:  (optional) function_handle specifying the density as a function of depth on the domain [-Lz 0]
            % - Parameter stratification:  (optional) function_handle specifying the stratification as a function of depth on the domain [-Lz 0]
            % - Parameter latitude: (optional) latitude of the domain (default is 33 degrees north)
            % - Parameter rho0: (optional) density at the surface z=0 (default is 1025 kg/m^3)
            % - Returns wvt: a new WVTransformHydrostatic instance
            arguments
                Lxyz (1,3) double {mustBePositive}
                Nxyz (1,3) double {mustBePositive}
                options.rho function_handle = @isempty
                options.N2 function_handle = @isempty
                options.dLnN2func function_handle = @isempty
                options.latitude (1,1) double = 33
                options.rho0 (1,1) double {mustBePositive} = 1025
                options.shouldAntialias logical = true
                options.jAliasingFraction double {mustBePositive(options.jAliasingFraction),mustBeLessThanOrEqual(options.jAliasingFraction,1)} = 2/3

                % ALL of these must be set for direct initialization to
                % avoid actually computing the modes.
                options.dLnN2 (:,1) double
                options.PFinv
                options.QGinv
                options.PF
                options.QG
                options.h (:,1) double
                options.P (:,1) double
                options.Q (:,1) double
                options.z (:,1) double
            end

            optionArgs = namedargs2cell(options);
            self@WVTransformAbstractHydrostatic(wvt,optionArgs{:});

            self.initializeStratifiedFlow();
            self.initializeGeostrophicComponent();
            self.initializeMeanDensityAnomalyComponent();
            self.initializeInternalGravityWaveComponent();
            self.initializeInertialOscillationComponent();
        end

        function wvtX2 = waveVortexTransformWithResolution(self,m)
            if ~isempty(self.dLnN2Function)
                wvtX2 = WVTransformHydrostatic([self.Lx self.Ly self.Lz],m, self.rhoFunction,latitude=self.latitude,rho0=self.rho0, N2func=self.N2Function, dLnN2func=self.dLnN2Function);
            else
                wvtX2 = WVTransformHydrostatic([self.Lx self.Ly self.Lz],m,latitude=self.latitude,rho0=self.rho0, N2=self.N2Function);
            end

            wvtX2.t0 = self.t0;
            wvtX2.t = self.t;
            [wvtX2.Ap,wvtX2.Am,wvtX2.A0] = self.spectralVariableWithResolution(wvtX2,self.Ap,self.Am,self.A0);
        end

        function h_0 = get.h_0(self)
            h_0 = self.h;
        end

        function h_pm = get.h_pm(self)
            h_pm = self.h;
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Transformations TO the spatial domain
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function u_bar = transformFromSpatialDomainWithFio(self, u)
            u_bar = (self.PF0*u)./self.P0;
        end

        function u_bar = transformFromSpatialDomainWithFg(self, u)
            u_bar = (self.PF0*u)./self.P0;
        end

        function w_bar = transformFromSpatialDomainWithGg(self, w)
            w_bar = (self.QG0*w)./self.Q0;
        end

        function w_bar = transformWithG_wg(~, w_bar )
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Needed to add and remove internal waves from the model
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function ratio = maxFw(self,kMode,lMode,j)
            arguments
                self WVTransform {mustBeNonempty}
                kMode (:,1) double
                lMode (:,1) double
                j (:,1) double
            end
            ratio = self.P0(j+1);
        end

        function ratio = maxFg(self,kMode,lMode,j)
            arguments
                self WVTransform {mustBeNonempty}
                kMode (:,1) double
                lMode (:,1) double
                j (:,1) double
            end
            ratio = self.P0(j+1);
        end

        [ncfile,matFilePath] = writeToFile(wvt,path,variables,options)

        function flag = isequal(self,other)
            arguments
                self WVTransform
                other WVTransform
            end
            flag = isequal@WVTransformAbstractHydrostatic(self,other);
            flag = flag & isequal(self.Ap, other.Ap);
            flag = flag & isequal(self.Am, other.Am);
            flag = flag & isequal(self.A0, other.A0);
        end
    end

    methods (Static)
        wvt = waveVortexTransformFromFile(path,options)
    end

end



