classdef ModelDiagnostics < handle
    % ModelDiagnostics
    properties (Access = public)
        wavemodel
        wintersmodel
    end
    
    methods
        function self = ModelDiagnostics(varargin)
            if isa(varargin{1}, 'WintersMode')
                self.wintersmodel = varargin{1};
                self.wavemodel = self.wintersmodel.wavemodel;
            elseif isa(varargin{1}, 'InternalWaveModel')
                self.wavemodel = varargin{1};
            end
        end
         
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         pv = LinearPV(self, u, v 
        
    end
end