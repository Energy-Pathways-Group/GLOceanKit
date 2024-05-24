classdef WVInertialOscillationMethods < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here

    properties (Abstract,GetAccess=public, SetAccess=public)
        Ap,Am,A0
    end
    properties (Abstract,GetAccess=public, SetAccess=protected)
        z
        UAp,VAp
        UAm,VAm
    end
    methods (Abstract)
        addPrimaryFlowComponent(self,primaryFlowComponent)
        u_bar = transformFromSpatialDomainWithFio(self,u)
    end

    methods
        % function WVInertialOscilationMethods()
        %     self.addPrimaryFlowComponent(flowComponent);
        % end
        function addInertialMotions(self,u,v)
            % add inertial motions to existing inertial motions
            %
            % The amplitudes of the inertial part of the flow will be added
            % to the existing inertial part of the flow.
            %
            % ```matlab
            % U_io = 0.2;
            % Ld = wvt.Lz/5;
            % u_NIO = @(z) U_io*exp((z/Ld));
            % v_NIO = @(z) zeros(size(z));
            %
            % wvt.addInertialMotions(u_NIO,v_NIO);
            % ```
            %
            % It is important to note that because the WVTransform
            % de-aliases by default, you will not likely get exactly the
            % same function out that you put in. The high-modes are
            % removed.
            %
            % The new inertial motions are added to the existing inertial motions
            % - Topic: Initial conditions — Inertial Oscillations
            % - Declaration: addInertialMotions(self,u,v)
            % - Parameter u: function handle that takes a single argument, u(Z)
            % - Parameter v: function handle that takes a single argument, v(Z)
            self.Ap(:,1) = self.Ap(:,1) + self.transformFromSpatialDomainWithFio(u(self.z) - sqrt(-1)*v(self.z))/2;
            self.Am(:,1) = conj(self.Ap(:,1));
        end

        function initWithInertialMotions(self,u,v)
            % initialize with inertial motions
            %
            % Clears variables Ap,Am,A0 and then sets inertial motions.
            % 
            % ```matlab
            % U_io = 0.2;
            % Ld = wvt.Lz/5;
            % u_NIO = @(z) U_io*exp((z/Ld));
            % v_NIO = @(z) zeros(size(z));
            %
            % wvt.initWithInertialMotions(u_NIO,v_NIO);
            % ```
            %
            % It is important to note that because the WVTransform
            % de-aliases by default, you will not likely get exactly the
            % same function out that you put in. The high-modes are
            % removed.
            %
            % - Topic: Initial conditions — Inertial Oscillations
            % - Declaration: initWithInertialMotions(self,u,v)
            % - Parameter u: function handle that takes a single argument, u(Z)
            % - Parameter v: function handle that takes a single argument, v(Z)
            self.Ap = zeros(size(self.Ap));
            self.Am = zeros(size(self.Am));
            self.A0 = zeros(size(self.A0));
            self.setInertialMotions(u,v);
        end

        function setInertialMotions(self,u,v)
            % set inertial motions
            %
            % % Overwrites existing inertial motions with the new values.
            % Other components of the flow will remain unaffected.
            %
            % ```matlab
            % U_io = 0.2;
            % Ld = wvt.Lz/5;
            % u_NIO = @(z) U_io*exp((z/Ld));
            % v_NIO = @(z) zeros(size(z));
            %
            % wvt.setInertialMotions(u_NIO,v_NIO);
            % ```
            %
            % It is important to note that because the WVTransform
            % de-aliases by default, you will not likely get exactly the
            % same function out that you put in. The high-modes are
            % removed.
            %
            
            % - Topic: Initial conditions — Inertial Oscillations
            % - Declaration: setInertialMotions(self,u,v)
            % - Parameter u: function handle that takes a single argument, u(Z)
            % - Parameter v: function handle that takes a single argument, v(Z)
            self.Ap(:,1) = self.transformFromSpatialDomainWithFio(u(self.z) - sqrt(-1)*v(self.z))/2;
            self.Am(:,1) = conj(self.Ap(:,1));
        end

        function removeAllInertialMotions(self)
            % remove all inertial motions
            %
            % All inertial motions are removed. Other components of the flow will remain unaffected.
            %
            % - Topic: Initial conditions — Inertial Oscillations
            % - Declaration: removeAllInertialMotions()
            self.Ap(:,1) = 0;
            self.Am(:,1) = 0;
        end
    end


end