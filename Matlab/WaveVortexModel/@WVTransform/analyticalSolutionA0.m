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
    options.shouldAssumeConstantN (1,1) logical {mustBeMember(options.shouldAssumeConstantN,[0 1])} = 1
end
if isfield(options,'A0')
    [kIndex,lIndex,jIndex] = self.geostrophicCoefficientsFromGeostrophicModes(kMode, lMode, jMode, 0, 0);
    A0Amp = options.A0;
else
    [kIndex,lIndex,jIndex] = self.geostrophicCoefficientsFromGeostrophicModes(kMode, lMode, jMode, 0, 0);
    A0Amp = 2*self.A0(kIndex,lIndex,jIndex);
end

A = abs(A0Amp);
phi = angle(A0Amp);

m = self.j(jIndex)*pi/self.Lz;
k = self.k(kIndex);
l = self.l(lIndex);
h = self.N0^2/(self.g*m^2);
sign = -2*(mod(jMode,2) == 1)+1;
norm = sign*sqrt(2*self.g/self.Lz)/self.N0;

G = @(z) norm*sin(m*(z+self.Lz));
F = @(z) norm*h*m*cos(m*(z+self.Lz));

theta = @(x,y,t) k*x + l*y + phi;
u = @(x,y,z,t) A*(self.g*l/self.f)*sin( theta(x,y,t) ).*F(z);
v = @(x,y,z,t) -A*(self.g*k/self.f)*sin( theta(x,y,t) ).*F(z);
w = @(x,y,z,t) zeros(self.Nx,self.Ny,self.Nz);
eta = @(x,y,z,t) A*cos( theta(x,y,t) ).*G(z);
p = @(x,y,z,t) A*self.rho0*self.g*cos( theta(x,y,t) ).*F(z);

end