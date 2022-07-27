classdef WVFlowConstituents
    properties
        bitmask = 0
    end
    properties (Constant)
        % binary 'literals' are only valid starting in Matlab 2019b
        none                = 0b00000000

        inertialOscillation = 0b00000001
        surfaceGravityWave  = 0b00000010
        internalGravityWave = 0b00000100
        surfaceGeostrophic  = 0b00001000
        internalGeostrophic = 0b00010000
        meanDensityAnomaly  = 0b00100000
        
        inertial            = 0b00000001
        wave                = 0b00000110
        geostrophic         = 0b00011000

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

