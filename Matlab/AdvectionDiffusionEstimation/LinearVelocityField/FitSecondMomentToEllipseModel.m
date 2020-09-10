function [parameters,error] = FitSecondMomentToEllipseModel( Mxx, Myy, Mxy, t, model)

if strcmp(model,'strain-diffusive')
    kappaScale = 1.0;
    sScale = 2e-6;
    model_strain_diffusive = @(a) MomentTensorModelError( Mxx, Myy, Mxy, t, 0, a(1), a(2), a(3), sScale, kappaScale, 'Ellipse-Overlap' );
    options = optimset('TolFun',1e-3,'TolX',1e-2, 'MaxFunEvals', 5000);
    [a, error,~,~] = fminsearch( model_strain_diffusive, [log(3.5e-6/sScale), -33*pi/180, log(0.2/kappaScale)], options );
    parameters.kappa = kappaScale*exp(a(3));
    parameters.sigma = sScale*exp(a(1));
    parameters.theta = a(2);
end


end

