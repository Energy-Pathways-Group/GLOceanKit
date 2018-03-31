N0 = 5.2e-3;
f0 = 7.9431e-05;
g = 9.81;
Lz = 5000;
k = 0.000001;

if k*k > (N0*N0 - f0*f0)/(g*Lz)
    fprintf('Using hyperbolic solution\n');
    f = @(q) Lz*(N0*N0 - f0*f0) - (1./q).*(g*(k*k*Lz*Lz-q.*q)).*tanh(q);
    h_f = @(m) (N0*N0 - f0*f0)./(g*(k*k - m*m ));
elseif k*k < (N0*N0 - f0*f0)/(g*Lz)
    fprintf('Using trig solution\n');
    f = @(q) Lz*(N0*N0 - f0*f0) - (1./q).*(g*(k*k*Lz*Lz+q.*q)).*tan(q);
    h_f = @(m) (N0*N0 - f0*f0)./(g*(k*k + m*m ));
else
    fprintf('lambda=0 case');
    m = 0;
    h = Lz;
end

m = fzero(f, k*Lz)/Lz;
h = h_f(m);
fprintf('(m,h) = (%.2g, %.2g)\n', m, h);

% f = @(q) Lz*(N0*N0 - f0*f0) - (1./q).*(g*(k*k*Lz*Lz-q.*q)).*tanh(q);
% w = @(q) Lz*(N0*N0 - f0*f0) - (1./q).*(g*(k*k*Lz*Lz+q.*q)).*tan(q);
% q = linspace((1-0.01)*k*Lz,(1+0.01)*k*Lz,100);
% figure, plot(q,f(q)), hold on, plot(q,w(q))