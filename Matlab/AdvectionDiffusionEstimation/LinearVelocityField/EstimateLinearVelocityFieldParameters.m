function [parameters,B] = EstimateLinearVelocityFieldParameters( x, y, t, parametersToEstimate, dof)
% EstimateLinearVelocityFieldParameters This function estimates linear
% velocity field parameters strain, vorticity and divergence from a cluster
% of Lagrangian particles.
%
% Required inputs are,
%   x ? [nT nDrifters] x position in meters of the nDrifters
%   y - [nT nDrifters] y position in meters of the nDrifters
%   t - [nT 1] times in seconds of the observation times
%   parametersToEstimate - array of ModelParameter objects.
%   dof (optional) ? positive integer indicating number of
%   degrees-of-freedom to allow in time for each parameter.
%
% Outputs are,
%   parameters - struct containing,
%       [u0,v0,u1,v1,sigma_n,sigma_s,zeta,delta] ? these will be [nT 1]
%       unless dof=1, in which case they will be scalar values.
%   B - array [nT nSplines] containing the B-splines used for the fit.
%
% Note that this implements equations 18-20 in Oscroft, Sykulski & Early. A
% proper implementation would remove the artificial requirement that
% drifter observations occur at the same time.
%
% This function calls EstimateLinearVelocityFieldParametersBSplines after
% generating B-splines following the description in the manuscript. Note
% that setting dof=1 is equivalent to constant parameters.

if ~exist('dof','var') || isempty(dof) || dof == 1
    K = 1;
    t_knot = [t(1) t(end)];
else
    K = min([dof 4],[]); % Only go as high as a cubic spline.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % This treats the center of each interval as a data point
    t_data = linspace(t(1),t(end),dof+1).';
    t_knot = InterpolatingSpline.KnotPointsForPoints((t_data(1:end-1)+t_data(2:end))/2,K);
    t_knot(1:K) = t_data(1);
    t_knot((end-K+1):end) = t_data(end);
end
B = BSpline.Spline(t,t_knot,K,0);

parameters = EstimateLinearVelocityFieldParametersBSplines(x, y, t, parametersToEstimate, B);

end

