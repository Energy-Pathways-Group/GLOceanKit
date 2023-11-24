function [u,v,w,eta,p] = analyticalSolutionApm(self, kMode, lMode, jMode, omegasign, options)
% return the analytical solution of the wave mode
%
% Returns function handles of the form u=@(x,y,z,t)
%
% - Topic: Analytical solutions
% - Declaration: [k,l] = analyticalSolutionApm(self)
% - Parameter kMode: integer index, (k0 > -Nx/2 && k0 < Nx/2)
% - Parameter lMode: integer index, (l0 > -Ny/2 && l0 < Ny/2)
% - Parameter jMode: integer index, (j0 >= 1 && j0 <= nModes), unless k=l=j=0
% - Parameter A0: (optional) amplitude in m. Default will use the current value.
% - Parameter phi: (optional) phase in radians, (0 <= phi <= 2*pi)
% - Parameter u: fluid velocity u (m/s)
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
    omegasign (1,1) double {mustBeMember(omegasign,[-1 1])}
    options.Apm (:,1) double
    options.phi (:,1) double = 0
    options.u (:,1) double
    options.shouldAssumeConstantN (1,1) logical {mustBeMember(options.shouldAssumeConstantN,[0 1])} = 1
end
if isfield(options,'u')
    [kIndex,lIndex,jIndex,ApAmp,AmAmp] = waveCoefficientsFromWaveModes(self, kMode, lMode, jMode, options.phi, options.u, omegasign);
    if ApAmp > 0
        Apm = ApAmp;
        omegasign = 1;
    else
        Apm = AmAmp;
        omegasign = -1;
    end 
elseif isfield(options,'Apm')
    [kIndex,lIndex,jIndex,ApAmp,AmAmp] = waveCoefficientsFromWaveModes(self, kMode, lMode, jMode, options.phi, 1, omegasign);
    if ApAmp > 0
        Apm = ApAmp;
        omegasign = 1;
    else
        Apm = AmAmp;
        omegasign = -1;
    end 
else
    [kIndex,lIndex,jIndex,ApAmp,AmAmp] = waveCoefficientsFromWaveModes(self, kMode, lMode, jMode, options.phi, 1, omegasign);
    if ApAmp > 0
        Apm = 2*self.Ap(kIndex,lIndex,jIndex);
        omegasign = 1;
    else
        Apm = 2*self.Am(kIndex,lIndex,jIndex);
        omegasign = -1;
    end 
end

A = abs(Apm);
phi = angle(Apm);

m = self.j(jIndex)*pi/self.Lz;
k = self.k(kIndex);
l = self.l(lIndex);
modesign = -2*(mod(jMode,2) == 1)+1;
if self.isHydrostatic
    h = self.N0^2/(self.g*m^2);
    norm = modesign*sqrt(2*self.g/self.Lz)/self.N0;
else
    h = (self.N0^2-self.f^2)/(k^2 + l^2 + m^2)/self.g;
    norm = modesign*sqrt(2*self.g/((self.N0^2 -self.f^2)*self.Lz));
end

G = @(z) norm*sin(m*(z+self.Lz));
F = @(z) norm*h*m*cos(m*(z+self.Lz));

alpha=atan2(l,k);
K = sqrt( k^2 + l^2);
omega = omegasign*sqrt( self.g*h*(k^2+l^2) + self.f^2 );
f0OverOmega = self.f/omega;
kOverOmega = K/omega;

theta = @(x,y,t) k*x + l*y + omega*t + phi;
u = @(x,y,z,t) A*(cos(alpha)*cos( theta(x,y,t) ) + f0OverOmega*sin(alpha)*sin( theta(x,y,t) )).*F(z);
v = @(x,y,z,t) A*(sin(alpha)*cos( theta(x,y,t) ) - f0OverOmega*cos(alpha)*sin( theta(x,y,t) )).*F(z);
w = @(x,y,z,t) A*K*h*sin( theta(x,y,t) ).*G(z);
eta = @(x,y,z,t) -A*h*kOverOmega * cos( theta(x,y,t)  ).*G(z);
p = @(x,y,z,t) -self.rho0*self.g*A*h*kOverOmega * cos( theta(x,y,t) ).*F(z);
end