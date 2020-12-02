function [parameters,B] = EstimateLinearVelocityFieldParameters( x, y, t, parametersToEstimate, dof)

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

