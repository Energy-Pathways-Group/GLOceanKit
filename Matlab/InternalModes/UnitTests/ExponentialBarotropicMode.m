z = linspace(-5000, 0, 1000)';
im = InternalModesExponentialStratification([5.2e-3 1300], [-5000 0], z, 33);
im.upperBoundary = UpperBoundary.freeSurface;

g = 9.81;
kArray = 10.^linspace(-7,1,9);
kArray = 1e1;

% Two relevant scales.
sqrt(g*im.Lz);
im.N0*im.b;

for k = kArray
    if k > sqrt(im.N0^2 - im.f0^2)/sqrt(g*im.Lz)
       c = sqrt(g*im.Lz/((2*k*k*(g*im.Lz)/(im.N0^2 - im.f0^2))-1));
       c = sqrt(g*tanh(k*im.Lz)/k); % can't find roots past N0*b, so just use this approximation.
    else
       c = sqrt(g*im.Lz);
    end
    h_guess = c*c/g
    
    x = im.b*im.N0/c;
    
    epsilon = im.f0/im.N0;
    lambda = k*im.b;
    nu = @(x) sqrt( epsilon^2 * x.^2 + lambda^2 );
    s = @(x) x;
    x_nu = lambda/sqrt(5*5*exp(-2*im.Lz/im.b) - epsilon*epsilon);
    
    r = im.FindRootsInRange(nu,s,[0.95*x 1.05*x],x_nu);
    fprintf('Found %d roots.\n',length(r));
    if length(r)>0
        h = (im.b*im.N0./r(1)).^2/g
    end
end


% omega = 2*im.N0;
% nu = @(x) omega*x/im.N0;
% s = @(x) im.N0*x/im.N0;
% 
% c = sqrt(g*im.Lz);
% x = im.b*im.N0/c;
% x_nu = Inf;
% 
% r = im.FindRootsInRange(nu, s, [0.95*x 1.05*x], x_nu);
% fprintf('Found %d roots.\n',length(r));
% if length(r)>0
%     h = (im.b*im.N0./r(1)).^2/g
% end