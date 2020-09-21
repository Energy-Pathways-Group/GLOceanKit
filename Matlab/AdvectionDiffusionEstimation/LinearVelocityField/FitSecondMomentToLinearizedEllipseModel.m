function [parameters,error] = FitSecondMomentToLinearizedEllipseModel(Mxx, Myy, Mxy, t, model)
if strcmp(model,'diffusive')

    
elseif strcmp(model,'strain-diffusive')
    dMxx = Mxx(end)-Mxx(1);
    dMyy = Myy(end)-Myy(1);
    dMxy = Mxy(end)-Mxy(1);
    theta = atan2(2*dMxy,dMxx-dMyy)/2;
    
    rotate_xx = @(Mxx,Myy,Mxy) cos(theta)*cos(theta)*Mxx + sin(theta)*sin(theta)*Myy + 2*cos(theta)*sin(theta)*Mxy;
    rotate_yy = @(Mxx,Myy,Mxy) sin(theta)*sin(theta)*Mxx + cos(theta)*cos(theta)*Myy - 2*cos(theta)*sin(theta)*Mxy;
    
    Mxx0_rotated = rotate_xx(Mxx(1),Myy(1),Mxy(1));
    Myy0_rotated = rotate_yy(Mxx(1),Myy(1),Mxy(1));
    MxxT_rotated = rotate_xx(Mxx(end),Myy(end),Mxy(end));
    MyyT_rotated = rotate_yy(Mxx(end),Myy(end),Mxy(end));
    
    T = t(end)-t(1);

    sigma = log( (MxxT_rotated+Myy0_rotated)/(Mxx0_rotated+MyyT_rotated) )/T;
    kappa = (sigma/2)*( MxxT_rotated*MyyT_rotated - Mxx0_rotated*Myy0_rotated)/( MxxT_rotated-Mxx0_rotated - MyyT_rotated + Myy0_rotated);
    
%     parameters.kappa = (1/(2*T)) * ( MxxT_rotated*Myy0_rotated +Mxx0_rotated*MyyT_rotated - 2*Mxx0_rotated*Myy0_rotated)/(Mxx0_rotated + Myy0_rotated);
%     parameters.sigma = (1/T) * ( MxxT_rotated - Mxx0_rotated + Myy0_rotated - MyyT_rotated)/(Mxx0_rotated + Myy0_rotated);
%     parameters.sigma = (1/T) * ( MxxT_rotated - Mxx0_rotated + Myy0_rotated - MyyT_rotated)/(Mxx0_rotated + Myy0_rotated);
    parameters.sigma = sigma;
    parameters.kappa = kappa;
    parameters.theta = theta;
end
error = 0;
parameters.zeta = 0;
end