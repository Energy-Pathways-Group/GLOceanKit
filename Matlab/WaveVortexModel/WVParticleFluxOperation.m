classdef WVParticleFluxOperation < WVOperation
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        isXYOnly
    end

    methods
        function self = WVParticleFluxOperation(name,f,options)
            arguments
                name char {mustBeNonempty}
                f function_handle
                options.isXYOnly double {mustBeMember(options.isXYOnly,[0 1])} = 0 
            end
            outputVariables(1) = WVVariableAnnotation('u_p',{'x','y','z'},'m/s', 'velocity x-direction at particle positions');
            outputVariables(2) = WVVariableAnnotation('v_p',{'x','y','z'},'m/s', 'velocity y-direction at particle positions');
            if ~options.isXYOnly
                outputVariables(3) = WVVariableAnnotation('w_p',{'x','y','z'},'m/s', 'velocity z-direction at particle positions');
            end

            self@WVOperation(name,outputVariables,f);
            self.isXYOnly = options.isXYOnly;
        end

        function varargout = Compute(self,wvt,x,y,z)
            varargout = cell(1,self.nVarOut);
            [varargout{:}] = self.f(wvt,x,y,z);
        end
    end
end