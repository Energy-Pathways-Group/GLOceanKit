classdef ParticleFluxOperation < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name % ?
        f = []
        xyOnly
    end

    methods
        function self = ParticleFluxOperation(name,f,xyOnly)
            arguments
                name char {mustBeNonempty}
                f {mustBeNonempty}
                xyOnly double {mustBeMember(xyOnly,[0 1])} = 0 
            end
            self.name = name;
            self.f = f;
            self.xyOnly = xyOnly;
        end

        function varargout = Compute(self,wvt,x,y,z)
            varargout = cell(1,length(self.outputVariables));
            [varargout{:}] = self.f(wvt,x,y,z);
        end
    end
end