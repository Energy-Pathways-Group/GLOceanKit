function [parameters,error,B] = FitSecondMomentToTimeVaryingEllipseModel( Mxx, Myy, Mxy, t, model,dof,K)

t_knot = InterpolatingSpline.KnotPointsForSplines(t,K,dof);

% B is matrix of size(B)=[length(t) nSplines]
B = BSpline.Spline(t,t_knot,K,0);
nSplines = size(B,2);

if strcmp(model,'strain-diffusive')
    kappaScale = 0.25;
    sScale = 2e-6;
    thetaScale = 36/pi;
    
    sigma0 = 4.0e-6; kappa0 = 0.5; theta0 = 0*pi/180;
    m_alpha = zeros(nSplines,1);
    m_s = B\(log(sigma0/sScale)*ones(size(t))); % y = B*m_s
    m_theta = B\(theta0/thetaScale*ones(size(t)));
    m_kappa = B\(log(kappa0/kappaScale)*ones(size(t)));
    a0 = cat(1,m_s,m_theta,m_kappa);
    
    model_strain_diffusive = @(a) MomentTensorModelErrorForBSplines( Mxx, Myy, Mxy, t, m_alpha, a(1:nSplines), a((nSplines+1):(2*nSplines)), a((2*nSplines+1):(3*nSplines)), sScale, kappaScale, thetaScale, 'Ellipse-Overlap', B );
    options = optimset('TolFun',1e-3,'TolX',1e-4, 'MaxFunEvals', 5000,'Display','iter');
    [a, error,~,~] = fminsearch( model_strain_diffusive, a0, options );
    
    kappa = (kappaScale*exp(B*a((2*nSplines+1):(3*nSplines))));
    sigma = (sScale*exp(B*a(1:nSplines)));
    theta = thetaScale*B*a((nSplines+1):(2*nSplines));
    
    figure, subplot(1,3,1), plot(kappa), subplot(1,3,2), plot(sigma), subplot(1,3,3), plot(theta)
    
    parameters.kappa = kappa;
    parameters.sigma = sigma;
    parameters.theta = theta;
end


end

