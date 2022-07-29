classdef WVTracerFluxOperation < WVOperation
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name % ?
        isXYOnly
    end

    methods
        function self = WVTracerFluxOperation(name,options)
            arguments
                name char {mustBeNonempty} = 'TracerFluxNoDamping'
                options.isXYOnly double {mustBeMember(options.isXYOnly,[0 1])} = 0 
            end
            self.name = name;
            self.isXYOnly = options.isXYOnly;
        end

        function varargout = compute(self,wvt,phi)
            varargout = cell(1,1);
            if self.isXYOnly
                varargout{1} = -wvt.u .* wvt.diffX(phi) - wvt.v .* wvt.diffY(phi);
            else
                varargout{1} = -wvt.u .* wvt.diffX(phi) - wvt.v .* wvt.diffY(phi) - wvt.w .* wvt.diffZF(phi);
            end
        end
    end
end