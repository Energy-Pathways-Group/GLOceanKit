classdef ParticleFluxOperation < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name % ?
        f = []
        xyOnly
        nVarOut
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
            if self.xyOnly == 1
                self.nVarOut = 2;
            else
                self.nVarOut = 3;
            end
        end

        function varargout = Compute(self,wvt,x,y,z)
            varargout = cell(1,self.nVarOut);
            [varargout{:}] = self.f(wvt,x,y,z);
        end
    end
end