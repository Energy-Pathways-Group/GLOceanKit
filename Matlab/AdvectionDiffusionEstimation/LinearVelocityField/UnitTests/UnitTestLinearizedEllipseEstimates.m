sigma = 2e-6;
kappa = 1;
theta = 0*45*pi/180;

Mxx0 = 1e5;
Myy0 = 1e5;
Mxy0 = 0;

t=linspace(0,10*86400,30).';

[Mxx, Myy, Mxy] = MomentTensorEvolutionInStrainVorticityField(Mxx0, Myy0, Mxy0, t, 0, sigma, theta, kappa);

parameters = FitSecondMomentToLinearizedEllipseModel(Mxx, Myy, Mxy, t, 'strain-diffusive');

sigmaEst = log( (Mxx+Myy0)./(Myy+Mxx0) )./t;
kappaEst = (sigma/2).*(Mxx.*Myy-Mxx0*Myy0)./(Mxx-Mxx0-Myy+Myy0);
