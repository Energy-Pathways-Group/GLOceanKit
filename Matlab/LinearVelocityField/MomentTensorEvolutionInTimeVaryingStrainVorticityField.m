% Mxx0, Myy0, Mxy0 are the initial moments, in units of meters^2.
% t is the time vector, in seconds, of size(t)=[n 1]
% zeta is the vorticity, 1/seconds, of size(zeta) = [n 1]
% sigma is the strain, 1/seconds, of size(zeta) = [n 1]
% theta is the angle of the principal strain axis, a counter-clockwise rotation, of size(zeta) = [n 1]
% kappa is the diffusivity, meters^2/second, of size(zeta) = [n 1]
% Mxx, Myy, Mxy will be returned
function [Mxx, Myy, Mxy] = MomentTensorEvolutionInTimeVaryingStrainVorticityField( Mxx0, Myy0, Mxy0, t, zeta, sigma, theta, kappa )

    Mxx = zeros(size(t));
    Myy = zeros(size(t));
    Mxy = zeros(size(t));

    Mxx(1) = Mxx0;
    Mxy(1) = Mxy0;
    Myy(1) = Myy0;

    for iTime=1:(length(t)-1)
        t_sub = t(iTime:(iTime+1))-t(iTime);
        [Mxx1, Myy1, Mxy1] = MomentTensorEvolutionInStrainVorticityField( Mxx(iTime), Myy(iTime), Mxy(iTime), t_sub, zeta(iTime), sigma(iTime), theta(iTime), kappa(iTime) );
        Mxx(iTime+1) = Mxx1(2);
        Mxy(iTime+1) = Mxy1(2);
        Myy(iTime+1) = Myy1(2);
    end
end