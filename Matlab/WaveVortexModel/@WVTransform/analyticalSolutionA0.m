function [u,v,w,eta,p] = analyticalSolutionA0(self, kMode, lMode, jMode, options)
% return the analytical solution of the geostrophic mode
%
% Returns function handles of the form u=@(x,y,z,t)
%
% - Topic: Analytical solutions
% - Declaration: [k,l] = analyticalSolutionA0(self)
% - Parameter kMode: integer index, (k0 > -Nx/2 && k0 < Nx/2)
% - Parameter lMode: integer index, (l0 > -Ny/2 && l0 < Ny/2)
% - Parameter jMode: integer index, (j0 >= 1 && j0 <= nModes), unless k=l=j=0
% - Parameter A0: (optional) amplitude in m. Default will use the current value.
% - Parameter phi: (optional) phase in radians, (0 <= phi <= 2*pi)
% - Parameter eta: (optional) isopycnal displacement (m)
% - Returns u: fluid velocity, u = @(x,y,z,t)
% - Returns v: fluid velocity, v = @(x,y,z,t)
% - Returns w: fluid velocity, w = @(x,y,z,t)
% - Returns eta: isopycnal displacement, eta = @(x,y,z,t)
% - Returns p: pressure, p = @(x,y,z,t)
arguments
    self WVTransform {mustBeNonempty}
    kMode (:,1) double
    lMode (:,1) double
    jMode (:,1) double
    options.A0 (:,1) double
    options.phi (:,1) double = 0
    options.u (:,1) double
    options.shouldAssumeConstantN (1,1) logical {mustBeMember(options.shouldAssumeConstantN,[0 1])} = 1
end
if isfield(options,'u')
    [kIndex,lIndex,jIndex,A0Amp] = self.geostrophicCoefficientsFromGeostrophicModes(kMode, lMode, jMode, options.phi, options.u);
elseif isfield(options,'A0')
    [kIndex,lIndex,jIndex] = self.geostrophicCoefficientsFromGeostrophicModes(kMode, lMode, jMode, options.phi, 0);
    A0Amp = options.A0;
else
    [kIndex,lIndex,jIndex] = self.geostrophicCoefficientsFromGeostrophicModes(kMode, lMode, jMode, options.phi, 0);
    A0Amp = wvt.A0(kIndex,lIndex,jIndex);
end

m = self.j(jIndex)*pi/self.Lz;
k = self.k(kIndex);
l = self.l(lIndex);
alpha=atan2(l,k);
K = sqrt( k^2 + l^2);
if j0 == 0
    omega = f;
else
    if self.isHydrostatic ==1
        omega = thesign*sqrt( (K*K*N0*N0 + m*m*f*f)/(m*m) );
    else
        omega = thesign*sqrt( (K*K*N0*N0 + m*m*f*f)/(K*K+m*m) );
    end
end

f0OverOmega = f/omega;
kOverOmega = K/omega;

 theta = @(x,y,t) kk*x + ll*y + omega*t + phi;
 u = @(x,y,z,t) U*(cos(alpha)*cos( theta(x,y,t) ) + f0OverOmega*sin(alpha)*sin( theta(x,y,t) )).*cos(m*z);
 v = @(x,y,z,t) U*(sin(alpha)*cos( theta(x,y,t) ) - f0OverOmega*cos(alpha)*sin( theta(x,y,t) )).*cos(m*z);
 w = @(x,y,z,t) (U*K/m) * sin( theta(x,y,t) ) .* sin(m*z);
 if j0 == 0
     zeta_unit = zeros(size(Z));
 else
     zeta_unit = -(U/m) * kOverOmega * cos( theta ) .* sin(m*Z);
 end
 rho_prime_unit = -(rho0/g)*N0*N0*(U*K/m/omega) * cos( theta ) .* sin(m*Z);
 rho_unit = rho0 -(N0*N0*rho0/g)*reshape(wavemodel.z,1,1,[]) + rho_prime_unit;

end