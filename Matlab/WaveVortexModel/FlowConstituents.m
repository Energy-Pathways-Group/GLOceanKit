classdef FlowConstituents < uint8
    methods
        function self = FlowConstituents(varargin)
            for k = 1:length(varargin)
                if strcmp(varargin{k}, 'inertialOscillation')
                    self = bitor(self,FlowConstituents.inertialOscillation);
                elseif strcmp(varargin{k}, 'surfaceGravityWave')
                    self = bitor(self,FlowConstituents.surfaceGravityWave);
                elseif strcmp(varargin{k}, 'internalGravityWave')
                    self = bitor(self,FlowConstituents.internalGravityWave);
                elseif strcmp(varargin{k}, 'surfaceGeostrophic')
                    self = bitor(self,FlowConstituents.surfaceGeostrophic);
                elseif strcmp(varargin{k}, 'internalGeostrophic')
                    self = bitor(self,FlowConstituents.internalGeostrophic);
                elseif strcmp(varargin{k}, 'meanDensityAnomaly')
                    self = bitor(self,FlowConstituents.meanDensityAnomaly);
                elseif strcmp(varargin{k}, 'inertial')
                    self = bitor(self,FlowConstituents.inertial);
                elseif strcmp(varargin{k}, 'wave')
                    self = bitor(self,FlowConstituents.wave);
                elseif strcmp(varargin{k}, 'geostrophic')
                    self = bitor(self,FlowConstituents.geostrophic);
                elseif strcmp(varargin{k}, 'all')
                    self = bitor(self,FlowConstituents.all);
                end
            end
        end
    end

    enumeration
        none                    (0)

        inertialOscillation     (0b00000001)
        surfaceGravityWave      (0b00000010)
        internalGravityWave     (0b00000100)
        surfaceGeostrophic      (0b00001000)
        internalGeostrophic     (0b00010000)
        meanDensityAnomaly      (0b00100000)
        
        inertial                (0b00000001)
        wave                    (0b00000110)
        geostrophic             (0b00011000)

        everything                     (0b00111111)
    end    
end

