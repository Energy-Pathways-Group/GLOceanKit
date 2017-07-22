
methods = cell(5,1);
methods{1} = 'finiteDifference';
methods{2} = 'wkbSpectral';
im = InternalModes('constant', 'wkbSpectral' , 64);


im.upperBoundary = 'free_surface';
% im.upperBoundary = 'rigid_lid';
im.normalization = 'max_u';

N0 = 5.2e-3;
f0 = 7.9431e-05;
g = 9.81;
Lz = 5000;
k_star = sqrt( (N0*N0 - f0*f0) / (g*Lz) );

im.ShowRelativeErrorAtWavenumber(0.1*k_star);
im.ShowRelativeErrorAtWavenumber(k_star);
im.ShowRelativeErrorAtWavenumber(10*k_star);

% im.ShowRelativeErrorAtFrequency(4*im.f0);

return







N0 = 5.2e-3;
f0 = 7.9431e-05;
g = 9.81;
Lz = 5000;
k = 0.01;

f_hyper = @(m) m*(N0*N0 - f0*f0)./(g*(k*k-m*m)) - tanh(m*Lz);
m(1) = fzero(f_hyper, max(k-1/Lz,1/Lz));
h2(1) = (N0*N0 - f0*f0)./(g*(k*k - m(1)*m(1) ));
fprintf('For hyperbolic: (m,h) = (%.2g, %.2g)\n', m(1), h2(1));

f_trig = @(m) m*(N0*N0 - self.f0*self.f0)./(g*(k*k+m*m)) - tan(m*Lz);
m(1) = fzero(f_trig, max(k-1/Lz,1/Lz));
h2(1) = (N0*N0 - f0*f0)./(g*(k*k + m(1)*m(1) ));

m = linspace(k - 0.1*k,k + 0.1*k,100);
figure
plot(m,f_hyper(m)), hold on
plot(m,f_trig(m))