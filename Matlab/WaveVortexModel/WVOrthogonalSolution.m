classdef WVOrthogonalSolution < handle
    %Analytical solution for one degree-of-freedom in the model
    %
    % Each degree-of-freedom in the model is associated with an analytical
    % solution to the equations of motion. This class provides a mapping
    % between the analytical solution and its numerical representation.
    %
    % The amplitude and phase are real valued. kMode, lMode, jMode are the
    % indices (k,l,j) used to represent the solution mathematically. Some
    % solution types might not have k or l (and thus be nil), and j may start
    % at 0, unlike Matlab indexing.
    %
    % (u,v,w,eta,p) are function handles which take arguments @(x,y,z,t) and
    % return real values.
    %
    % The numerical representation of that solution includes the location of
    % the solution and its expected amplitude in Ap/Am/A0 coefficient matrix.
    % Solutions may have their conjugate located in the same (or different)
    % matrix.
    %
    % The energyFactor and enstrophyFactor indicate the multiplicative factor
    % required to multiply the sum of the squared amplitudes by to get energy
    % and enstrophy.
    %
    % - Declaration: classdef WVAnalyticalSolution < handle
    properties
        kMode, lMode, jMode
        amplitude
        phase
        u,v,w,eta,rho_e,p,ssh,qgpv,psi
        Lx,Ly,Lz
        N2

        coefficientMatrix WVCoefficientMatrix
        coefficientMatrixIndex
        coefficientMatrixAmplitude = 0

        conjugateCoefficientMatrix WVCoefficientMatrix
        conjugateCoefficientMatrixIndex
        conjugateCoefficientMatrixAmplitude = 0

        energyFactor
        enstrophyFactor
    end

    methods
        function self = WVOrthogonalSolution(kMode, lMode, jMode, amplitude,phase,u,v,w,eta,rho_e,p,ssh,qgpv,options)
            arguments (Input)
                kMode (1,1) double
                lMode (1,1) double
                jMode (1,1) double
                amplitude (1,1) double
                phase (1,1) double
                u function_handle
                v function_handle
                w function_handle
                eta function_handle
                rho_e function_handle
                p function_handle
                ssh function_handle
                qgpv function_handle
                options.psi function_handle = @(x,y,z,t) 0*x;
                options.Lxyz (3,1) double
                options.N2 function_handle
            end
            self.amplitude = amplitude;
            self.phase = phase;
            self.kMode = kMode;
            self.lMode = lMode;
            self.jMode = jMode;
            self.u = u;
            self.v = v;
            self.w = w;
            self.eta = eta;
            self.rho_e = rho_e;
            self.p = p;
            self.ssh = ssh;
            self.qgpv = qgpv;
            self.psi = options.psi;
            self.Lx = options.Lxyz(1);
            self.Ly = options.Lxyz(2);
            self.Lz = options.Lxyz(3);
            self.N2 = options.N2;
        end

        function energy = depthIntegratedTotalEnergy(self,options)
            arguments (Input)
                self WVOrthogonalSolution {mustBeNonempty}
                options.isHydrostatic (1,1) logical {mustBeMember(options.isHydrostatic,[0 1])} = 0
            end
            arguments (Output)
                energy (1,1) double
            end
            if options.isHydrostatic==1
                energyDensity = @(x,y,z) self.u(x,y,z,0).*self.u(x,y,z,0) + self.v(x,y,z,0).*self.v(x,y,z,0) + self.N2(z).*self.eta(x,y,z,0).*self.eta(x,y,z,0);
            else
                energyDensity = @(x,y,z) self.u(x,y,z,0).*self.u(x,y,z,0) + self.v(x,y,z,0).*self.v(x,y,z,0) + self.w(x,y,z,0).*self.w(x,y,z,0) + self.N2(z).*self.eta(x,y,z,0).*self.eta(x,y,z,0);
            end
            energy = integral3(energyDensity,0,self.Lx,0,self.Ly,-self.Lz,0,AbsTol=1e-12)/self.Lx/self.Ly/2;
        end

        function enstrophy = depthIntegratedTotalEnstrophy(self)
            arguments (Input)
                self WVOrthogonalSolution {mustBeNonempty}
            end
            arguments (Output)
                enstrophy (1,1) double
            end
            enstrophy = integral3(@(x,y,z) self.qgpv(x,y,z,0).*self.qgpv(x,y,z,0),0,self.Lx,0,self.Ly,-self.Lz,0,AbsTol=1e-12)/self.Lx/self.Ly/2;
        end
    end

end