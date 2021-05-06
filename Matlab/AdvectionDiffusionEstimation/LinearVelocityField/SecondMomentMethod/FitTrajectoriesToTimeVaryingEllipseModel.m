function [parameters,error,B] = FitTrajectoriesToTimeVaryingEllipseModel( x, y, t, model, dof, K)

[~, ~, q, r] = CenterOfMass( x, y );
[Mxx, Myy, Mxy] = SecondMomentMatrix( q, r);

[parameters,error, B] = FitSecondMomentToTimeVaryingEllipseModel( Mxx, Myy, Mxy, t, model,dof,K);

end

