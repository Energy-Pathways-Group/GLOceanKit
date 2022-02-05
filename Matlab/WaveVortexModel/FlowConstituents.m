classdef FlowConstituents
    properties
        bitmask = 0
    end
    properties (Constant)
        none                = 0b00000001

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
        function self = FlowConstituents(varargin)
            bitmask = 0;
            for k = 1:length(varargin)
                if strcmp(varargin{k}, 'inertialOscillation')
                    bitmask = bitor(bitmask,FlowConstituents.inertialOscillation);
                elseif strcmp(varargin{k}, 'surfaceGravityWave')
                    bitmask = bitor(bitmask,FlowConstituents.surfaceGravityWave);
                elseif strcmp(varargin{k}, 'internalGravityWave')
                    bitmask = bitor(bitmask,FlowConstituents.internalGravityWave);
                elseif strcmp(varargin{k}, 'surfaceGeostrophic')
                    bitmask = bitor(bitmask,FlowConstituents.surfaceGeostrophic);
                elseif strcmp(varargin{k}, 'internalGeostrophic')
                    bitmask = bitor(bitmask,FlowConstituents.internalGeostrophic);
                elseif strcmp(varargin{k}, 'meanDensityAnomaly')
                    bitmask = bitor(bitmask,FlowConstituents.meanDensityAnomaly);
                elseif strcmp(varargin{k}, 'inertial')
                    bitmask = bitor(bitmask,FlowConstituents.inertial);
                elseif strcmp(varargin{k}, 'wave')
                    bitmask = bitor(bitmask,FlowConstituents.wave);
                elseif strcmp(varargin{k}, 'geostrophic')
                    bitmask = bitor(bitmask,FlowConstituents.geostrophic);
                elseif strcmp(varargin{k}, 'all')
                    bitmask = bitor(bitmask,FlowConstituents.all);
                end
            end
        end
        
        function bool = Contains(self,otherFlowConstituent)
            if isa(otherFlowConstituent,"numeric")
                bool = logical(bitor(self.bitmask,otherFlowConstituent));
            elseif isa(otherFlowConstituent,"FlowConstituents")
                bool = logical(bitor(self.bitmask,otherFlowConstituent.bitmask));
            end
        end

    end 
end

