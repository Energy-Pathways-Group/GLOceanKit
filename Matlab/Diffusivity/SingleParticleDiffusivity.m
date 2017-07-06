function [r2, kappa_a ] = SingleParticleDiffusivity( t, x, y, method )
%SingleParticleDiffusivity Compute the single particle diffusivity
% *relative to the center of mass* of the given trajectories.
%
% Note that the mean-square-distance, r2, can bias a sample. If you start a
% bunch of particles at the origin on a random walk, you will find a
% correlation with diffusivity and mean-square-distance because you've
% selectively chosen the particles that diffused the furthest.
nDrifters = size(x,2);
r2 = zeros(nDrifters,1);
kappa_a = zeros(nDrifters,1);
velocities = zeros(length(t)-1,nDrifters);
dt = t(2)-t(1);

[~, ~, xc, yc] = CenterOfMass( x, y );

stride = 1;
for iDrifter=1:stride:nDrifters
    q = xc(:,iDrifter);
    r = yc(:,iDrifter);
    
    % mean-squared-separation distance over the requested time-interval
    r2(iDrifter) = mean(q.^2 + r.^2,1);
    
    % Now remove the initial conditions.
    q = q-q(1);
    r = r-r(1);
    
    % squared-separation distance as a function of time.
    D2 = q.^2 + r.^2;
    
    if strcmp(method,'endpoint')
        kappa_a(iDrifter) = (D2(end) - D2(1))/(4*(t(end) - t(1)));
    elseif strcmp(method,'slope')
        [p,~,mu]=polyfit(t,D2,1);
        kappa_a(iDrifter) = (p(1)/mu(2))/4;
    elseif strcmp(method,'powspec')
        velocities(:,iDrifter) = diff(q+sqrt(-1)*r)/dt;
    elseif strcmp(method,'y0v')
        % This is not a complete measure of diffusivity---this is a
        % term that we hope goes to zero. Note that if you take the
        % mean of the relative velocity (in time), that will be
        % positive, if the particles are separating on average.
        kappa_a(iDrifter) = q(1)*(q(end)-q(end-1))/dt + r(1)*(r(end)-r(end-1))/dt;
    elseif strcmp(method,'yv')
        % This is identical to the 'endpoint' method.
        u = diff(q)/dt;
        v = diff(r)/dt;
        qc = q(2:end) - diff(q)/2;
        rc = r(2:end) - diff(r)/2;
        kappa_a(iDrifter) = (mean(u.*qc) + mean(v.*rc))/2;
    end
end

if strcmp(method,'powspec')
    averaging_bandwidth = 1;
    taper_bandwidth = 0;
    kappa_a = DiffusivityFromZeroFrequency(t(2)-t(1),velocities,averaging_bandwidth,taper_bandwidth)';
end

end

