classdef APEOperation < WVOperation
% Computes the nonlinear flux for a WVTransform

    properties (GetAccess=public, SetAccess=protected)
    
    end

    methods

        function self = APEOperation()
            outputVariables(1) = WVVariableAnnotation('ape',{'x','y','z'},'m^2 s^{-2}', 'available potential energy');
            self@WVOperation('ape',outputVariables,@disp);
        end

        function varargout = compute(self,wvt,varargin)
            eta_true = wvt.eta_true;
            ape = zeros(wvt.spatialMatrixSize);
            Z = wvt.Z;
            for i=1:numel(Z)
                z = Z(i);
                eta = eta_true(i);
                ape(i) = integral( @(xi) xi.*wvt.N2Function(z-xi),0,eta);
            end
            varargout = {ape};
        end
  
    end

end