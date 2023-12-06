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
    u,v,w,eta,p

    coefficientMatrix WVCoefficientMatrix
    coefficientMatrixIndex
    coefficientMatrixAmplitude

    conjugateCoefficientMatrix WVCoefficientMatrix
    conjugateCoefficientMatrixIndex
    conjugateCoefficientMatrixAmplitude

    energyFactor
    enstrophyFactor
end

methods
    function self = WVOrthogonalSolution(kMode, lMode, jMode, amplitude,phase,u,v,w,eta,p)
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
            p function_handle
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
        self.p = p;
    end
end

end