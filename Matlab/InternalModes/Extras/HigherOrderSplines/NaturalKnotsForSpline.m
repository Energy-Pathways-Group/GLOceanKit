%% NaturalKnotsForSpline
%
% Returns the natural knot points for splines of order K. It defaults to 1
% degree of freedom, but you can override this by setting DF to some
% nonnegative integer value.
function [t_knot] = NaturalKnotsForSpline( t, K, DF )

if nargin < 3
    DF = 1;
end

if (DF < 1 || mod(DF,1) ~= 0)
    disp('DF must be a non-negative integer');
    return;
end

t = [t(1); t(1+DF:DF:end-DF); t(end)];

if mod(K,2) == 1
    % Odd spline order, so knots go in between points.
    dt = diff(t);
    
    % This gives us N+1 knot points
    t_knot = [t(1); t(1:end-1)+dt/2; t(end)];
    
    % Now remove start and end knots
    for i=1:((K-1)/2)
        t_knot(2) = [];
        t_knot(end-1) = [];
    end
    
else
    t_knot = t;
    
    % Now remove start and end knots
    for i=1:((K-2)/2)
        t_knot(2) = [];
        t_knot(end-1) = [];
    end
    
end

end