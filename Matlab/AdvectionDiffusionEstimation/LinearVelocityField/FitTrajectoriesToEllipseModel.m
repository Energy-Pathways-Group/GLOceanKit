function [parameters,error] = FitTrajectoriesToEllipseModel( x, y, t, model)

[~, ~, q, r] = CenterOfMass( x, y );
[Mxx, Myy, Mxy] = SecondMomentMatrix( q, r);

[parameters,error] = FitSecondMomentToEllipseModel( Mxx, Myy, Mxy, t, model, 'fminsearch');

% [parameters1,error1] = FitSecondMomentToEllipseModel( Mxx, Myy, Mxy, t, model, 'fminsearch');
% [parameters2,error2] = FitSecondMomentToEllipseModel( Mxx, Myy, Mxy, t, model, 'angle-gridded');
% [parameters3,error3] = FitSecondMomentToEllipseModel( Mxx, Myy, Mxy, t, model, 'gridded');
% 
% r13 = 100*(error1/error3-1);
% r23 = 100*(error2/error3-1);
% 
% fprintf('Compared to the gridded, fminsearch and angle-gridded performed %.1f and %.1f percent.\n',r13,r23);
% parameters = parameters3;
% error = error3;


end

