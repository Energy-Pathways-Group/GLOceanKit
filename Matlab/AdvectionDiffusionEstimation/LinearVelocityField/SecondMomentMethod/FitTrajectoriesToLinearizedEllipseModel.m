function [parameters,error] = FitTrajectoriesToLinearizedEllipseModel( x, y, t, model, varargin)

[~, ~, q, r] = CenterOfMass( x, y );
[Mxx, Myy, Mxy] = SecondMomentMatrix( q, r);
[parameters,error] = FitSecondMomentToLinearizedEllipseModel( Mxx, Myy, Mxy, t, model, varargin{:});

end