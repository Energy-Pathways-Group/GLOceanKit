z = linspace(-5000, 0, 1000)';
im = InternalModesExponentialStratification([5.2e-3 1300], [-5000 0], z, 33);
im.upperBoundary = UpperBoundary.freeSurface;

g = 9.81;


% Two relevant scales.
% sqrt(g*im.Lz);
% im.N0*im.b;

f = @(k) im.b*im.N0./sqrt(g*tanh(k*im.Lz)./k);

kArray = 10.^linspace(-7,-1,9);
xFound = zeros(size(kArray));
for i = 1:length(kArray)
    k = kArray(i)
    x = f(k);
    
    epsilon = im.f0/im.N0;
    lambda = k*im.b;
    nu = @(x) sqrt( epsilon^2 * x.^2 + lambda^2 );
    s = @(x) x;
    x_nu = lambda/sqrt(5*5*exp(-2*im.Lz/im.b) - epsilon*epsilon);
    
    r = im.FindRootsInRange(nu,s,[0.95*x 1.05*x],x_nu);

    if ~isempty(r)
        xFound(i) = r(1);
    end
end


figure
plot(kArray,f(kArray)), hold on
scatter(kArray,xFound)


x_nu = Inf;
f = @(omega) max( im.b*im.N0*omega/g, im.b*im.N0/sqrt(g*im.Lz));

omegaArray = im.N0*linspace(0.5,30,100)';
xFound = zeros(size(omegaArray));

for i = 1:length(omegaArray)
    omega = omegaArray(i);
    x = f(omega);
    nu = @(x) omega*x/im.N0;
    s = @(x) im.N0*x/im.N0;
    r = im.FindRootsInRange(nu, s, [0.95*x 1.5*x], x_nu);
    if ~isempty(r)
        xFound(i) = r(1);
    end
end

figure
plot(omegaArray,f(omegaArray)), hold on
scatter(omegaArray,xFound)