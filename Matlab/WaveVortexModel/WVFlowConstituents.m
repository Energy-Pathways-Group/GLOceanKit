classdef WVFlowConstituents
    %Describes a set of flow constituents
    %
    % The WVFlowConstituents class is used to describe sets of flow
    % constituents. This class is used by methods that extract only part of
    % the flow.
    %
    % As an example, if you initialize a variable called `flowConstituent`
    % with
    % ```matlab
    % flowConstituent = WVFlowConstituents('geostrophic')
    % ```
    % or equivalently with,
    % ```matlab
    % flowConstituent = WVFlowConstituents.geostrophic
    % ```
    % then you could call
    % ```matlab
    % [u,v,w] = wvt.velocityOfFlowConstituent(flowConstituent)
    % ```
    % and that would return only the geostrophic portion of the flow.
    %
    % - Declaration: classdef WVFlowConstituents
    properties (Access=private)
        bitmask = 0
    end
    properties (Constant)
        % no constituents (no flow)
        % - Topic: Constituent groups
        none                = 0b00000000

        % inertial oscillations
        % - Topic: Primary constituents
        inertialOscillation = 0b00000001

        % surface gravity waves
        % - Topic: Primary constituents
        surfaceGravityWave  = 0b00000010

        % internal gravity waves
        % - Topic: Primary constituents
        internalGravityWave = 0b00000100

        % surface geostrophic motions
        % - Topic: Primary constituents
        surfaceGeostrophic  = 0b00001000

        % internal geostrophic motions
        % - Topic: Primary constituents
        internalGeostrophic = 0b00010000

        % mean density anomaly
        % - Topic: Primary constituents
        meanDensityAnomaly  = 0b00100000

        % inertial oscillations
        % - Topic: Constituent groups
        inertial            = 0b00000001

        % surface and internal gravity waves
        % - Topic: Constituent groups
        wave                = 0b00000110

        % surface and internal geostrophic motions
        % - Topic: Constituent groups
        geostrophic         = 0b00011000

        % all motions (all flow)
        % - Topic: Constituent groups
        all                 = 0b00111111
    end
    methods
        function self = WVFlowConstituents(varargin)
            self.bitmask = 0;
            for k = 1:length(varargin)
                if strcmp(varargin{k}, 'inertialOscillation')
                    self.bitmask = bitor(self.bitmask,WVFlowConstituents.inertialOscillation);
                elseif strcmp(varargin{k}, 'surfaceGravityWave')
                    self.bitmask = bitor(self.bitmask,WVFlowConstituents.surfaceGravityWave);
                elseif strcmp(varargin{k}, 'internalGravityWave')
                    self.bitmask = bitor(self.bitmask,WVFlowConstituents.internalGravityWave);
                elseif strcmp(varargin{k}, 'surfaceGeostrophic')
                    self.bitmask = bitor(self.bitmask,WVFlowConstituents.surfaceGeostrophic);
                elseif strcmp(varargin{k}, 'internalGeostrophic')
                    self.bitmask = bitor(self.bitmask,WVFlowConstituents.internalGeostrophic);
                elseif strcmp(varargin{k}, 'meanDensityAnomaly')
                    self.bitmask = bitor(self.bitmask,WVFlowConstituents.meanDensityAnomaly);
                elseif strcmp(varargin{k}, 'inertial')
                    self.bitmask = bitor(self.bitmask,WVFlowConstituents.inertial);
                elseif strcmp(varargin{k}, 'wave')
                    self.bitmask = bitor(self.bitmask,WVFlowConstituents.wave);
                elseif strcmp(varargin{k}, 'geostrophic')
                    self.bitmask = bitor(self.bitmask,WVFlowConstituents.geostrophic);
                elseif strcmp(varargin{k}, 'all')
                    self.bitmask = bitor(self.bitmask,WVFlowConstituents.all);
                end
            end
        end
        
        function bool = Contains(self,otherFlowConstituent)
            if isa(otherFlowConstituent,"numeric")
                bool = logical(bitand(self.bitmask,otherFlowConstituent));
            elseif isa(otherFlowConstituent,"WVFlowConstituents")
                bool = logical(bitand(self.bitmask,otherFlowConstituent.bitmask));
            end
        end

    end 
end

