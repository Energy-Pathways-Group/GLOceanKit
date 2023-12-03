classdef WVAnalyticalSolution < handle
%Analytical solution for one degree-of-freedom in the model
% 
% Each degree-of-freedom in the model is associated with an analytical
% solution to the equations of motion. This class provides a mapping
% between the analytical solution and its numerical representation.
%
% The amplitude and phase are real valued. kMode, lMode, jMode are the
% indices (k,l,j) used to represent the solution mathematically. Some
% solution types might not have k or l (and thus be nil), and j may start
% at 0, unlike Matlab index.
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
    amplitude
    phase
    kMode, lMode, jMode
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

end