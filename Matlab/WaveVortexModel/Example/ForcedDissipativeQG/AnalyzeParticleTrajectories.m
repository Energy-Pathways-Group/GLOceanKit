ncfile = NetCDFFile('ForcedDissipativeQG-particles-512.nc');
[x,y,eta,t] = ncfile.readVariables('drifter-x','drifter-y','drifter-eta','t');
x = x.';
y = y.';

rms_eta = sqrt(mean(eta.^2,2,'omitnan'));
etaAvg = mean(eta,1,'omitnan').';
figure
subplot(1,2,1)
histogram((rms_eta.^(-1/3))/86400)
ylabel('rms enstrophy cascade time scale per trajectory (days)')
subplot(1,2,2)
plot(t/86400,((-etaAvg).^(-1/3))/86400)
xlabel('t (days)')
ylabel('spatial mean enstrophy cascade time scale (days)')

enstrophyTimeScale = mean(((-etaAvg).^(-1/3)));


cx = x+sqrt(-1)*y;
dt = t(2)-t(1);
cv = vdiff( dt, cx, 1);

cv = cv(:,1:100:end);

[psi,lambda]=sleptap(size(cv,1),3);
[f,spp,snn,spn]=mspec(dt, cv,psi);

params = maternfit(dt,cv,1/enstrophyTimeScale);

% The integral time scale of the Matern is found with equation 43 in the
% Lilly, et al. NPG paper. Take S(0)/sigma^2 to get, 1/(\lambda * c_alpha)
% where c_alpha = beta(1/2, alpha-1/2)/(2*pi)
T_decorrelation = 1./(params.lambda .* beta(1/2,params.alpha-1/2));

i = 20;
[f,s] = maternspec(dt,1000,params.sigma(i),params.alpha(i),params.lambda(i));
std(cv(:,i))^2
pi*params.sigma(i)^2
sum(s)*(f(2)-f(1))

% matern = @(omega,A,lambda,alpha) A^2./(omega.^2 + lambda^2).^alpha;
% i = 1;
% f = @(omega) matern(omega,params.)

figure, scatter(T_decorrelation,(rms_eta(1:100:end)).^(-1/3))


