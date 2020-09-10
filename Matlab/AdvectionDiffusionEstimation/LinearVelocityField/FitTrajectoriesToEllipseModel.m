function [parameters,error] = FitTrajectoriesToEllipseModel( x, y, t, model)

[~, ~, q, r] = CenterOfMass( x, y );
[Mxx, Myy, Mxy] = SecondMomentMatrix( q, r);

[parameters,error] = FitSecondMomentToEllipseModel( Mxx, Myy, Mxy, t, model);

end

