classdef FlowConstituents
    properties
        bitmask = 0
    end
    properties (Constant)
        none                = 0

        inertialOscillation = 1
        surfaceGravityWave  = 2
        internalGravityWave = 4
        surfaceGeostrophic  = 8
        internalGeostrophic = 16
        meanDensityAnomaly  = 32
        
        inertial            = 1
        wave                = 2+4
        geostrophic         = 8+16

        all                 = 1+2+4+8+16+32


        % binary 'literals' are only valid starting in Matlab 2019b
%         none                = 0b00000000
% 
%         inertialOscillation = 0b00000001
%         surfaceGravityWave  = 0b00000010
%         internalGravityWave = 0b00000100
%         surfaceGeostrophic  = 0b00001000
%         internalGeostrophic = 0b00010000
%         meanDensityAnomaly  = 0b00100000
%         
%         inertial            = 0b00000001
%         wave                = 0b00000110
%         geostrophic         = 0b00011000
% 
%         all                 = 0b00111111
    end
    methods
        function self = FlowConstituents(varargin)
            self.bitmask = 0;
            for k = 1:length(varargin)
                if strcmp(varargin{k}, 'inertialOscillation')
                    self.bitmask = bitor(self.bitmask,FlowConstituents.inertialOscillation);
                elseif strcmp(varargin{k}, 'surfaceGravityWave')
                    self.bitmask = bitor(self.bitmask,FlowConstituents.surfaceGravityWave);
                elseif strcmp(varargin{k}, 'internalGravityWave')
                    self.bitmask = bitor(self.bitmask,FlowConstituents.internalGravityWave);
                elseif strcmp(varargin{k}, 'surfaceGeostrophic')
                    self.bitmask = bitor(self.bitmask,FlowConstituents.surfaceGeostrophic);
                elseif strcmp(varargin{k}, 'internalGeostrophic')
                    self.bitmask = bitor(self.bitmask,FlowConstituents.internalGeostrophic);
                elseif strcmp(varargin{k}, 'meanDensityAnomaly')
                    self.bitmask = bitor(self.bitmask,FlowConstituents.meanDensityAnomaly);
                elseif strcmp(varargin{k}, 'inertial')
                    self.bitmask = bitor(self.bitmask,FlowConstituents.inertial);
                elseif strcmp(varargin{k}, 'wave')
                    self.bitmask = bitor(self.bitmask,FlowConstituents.wave);
                elseif strcmp(varargin{k}, 'geostrophic')
                    self.bitmask = bitor(self.bitmask,FlowConstituents.geostrophic);
                elseif strcmp(varargin{k}, 'all')
                    self.bitmask = bitor(self.bitmask,FlowConstituents.all);
                end
            end
        end
        
        function bool = Contains(self,otherFlowConstituent)
            if isa(otherFlowConstituent,"numeric")
                bool = logical(bitand(self.bitmask,otherFlowConstituent));
            elseif isa(otherFlowConstituent,"FlowConstituents")
                bool = logical(bitand(self.bitmask,otherFlowConstituent.bitmask));
            end
        end

    end 
end

