classdef WVInertialOscillationMethods < handle
    %UNTITLED2 Summary of this class goes here
    %   Detailed explanation goes here


    methods
        function addInertialMotions(self,u,v)
            % add inertial motions to existing inertial motions
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
            % The new inertial motions are added to the existing inertial motions
            % - Topic: Initial conditions — Inertial Oscillations
            % - Declaration: addInertialMotions(self,u,v)
            % - Parameter u: function handle that takes a single argument, u(Z)
            % - Parameter v: function handle that takes a single argument, v(Z)
            self.Ap(:,1) = self.transformFromSpatialDomainWithFio(u(self.z) - sqrt(-1)*v(self.z))/2;
            self.Am(:,1) = conj(self.Ap(:,1));
        end

        function initWithInertialMotions(self,u,v)
            % initialize with inertial motions
            %
            % ```matlab
            % U_io = 0.2;
            % Ld = wvt.Lz/5;
            % u_NIO = @(z) U_io*exp(-(z/Ld));
            % v_NIO = @(z) zeros(size(z));
            %
            % wvt.initWithInertialMotions(u_NIO,v_NIO);
            % ```
            %
            % Clears variables Ap,Am,A0 and then sets inertial motions
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
            % ```matlab
            % U_io = 0.2;
            % Ld = wvt.Lz/5;
            % u_NIO = @(z) U_io*exp(-(z/Ld));
            % v_NIO = @(z) zeros(size(z));
            %
            % wvt.setInertialMotions(u_NIO,v_NIO);
            % ```
            %
            % Overwrites existing inertial motions with the new values
            % - Topic: Initial conditions — Inertial Oscillations
            % - Declaration: setInertialMotions(self,u,v)
            % - Parameter u: function handle that takes a single argument, u(Z)
            % - Parameter v: function handle that takes a single argument, v(Z)
            [~,~,Z] = self.xyzGrid;

            Ubar = self.transformFromSpatialDomainWithF( u(Z) );
            Vbar = self.transformFromSpatialDomainWithF( v(Z) );
            Ap_ = self.ApU.*Ubar + self.ApV.*Vbar;
            Am_ = self.AmU.*Ubar + self.AmV.*Vbar;
            self.Ap(1,1,:) = Ap_(1,1,:);
            self.Am(1,1,:) = Am_(1,1,:);
        end

        function removeAllInertialMotions(self)
            % remove all inertial motions
            %
            % All inertial motions are removed
            % - Topic: Initial conditions — Inertial Oscillations
            % - Declaration: removeAllInertialMotions()
            self.Ap(1,1,:) = 0;
            self.Am(1,1,:) = 0;
        end
    end


end