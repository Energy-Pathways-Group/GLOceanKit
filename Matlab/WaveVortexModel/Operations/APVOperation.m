classdef APVOperation < WVOperation
% Computes the nonlinear flux for a WVTransform

    properties (GetAccess=public, SetAccess=protected)
    
    end

    methods

        function self = APVOperation()
            outputVariables(1) = WVVariableAnnotation('apv',{'x','y','z'},'1/s', 'available potential vorticity');
            self@WVOperation('apv',outputVariables,@disp);
        end

        function varargout = compute(self,wvt,varargin)
            eta_true = wvt.eta_true;
            % zeta_x = wvt.diffY(wvt.w) - wvt.diffZF(wvt.v); % w_y - v_z
            % zeta_y = wvt.diffZF(wvt.u) - wvt.diffX(wvt.w);  % u_z - w_x
            % zeta_z = wvt.diffX(wvt.v) - wvt.diffY(wvt.u);  % v_x - u_y

            APV_x = wvt.zeta_x .* wvt.diffX(eta_true);
            APV_y = wvt.zeta_y .* wvt.diffY(eta_true);
            APV_z = (wvt.zeta_z + wvt.f) .* wvt.diffZG(eta_true);
            APV = wvt.zeta_z - APV_x - APV_y - APV_z;
            varargout = {APV};
        end
  
    end

end