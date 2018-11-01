function [r2, kappa_r, velocity_correlation ] = PairwiseRelativeDiffusivity( t, x, y, method )
%PairwiseRelativeDiffusivity Compute the relative diffusivity using all
% pairwise combinations of particles.
%
% Note that the mean-square-distance, r2, can bias a sample. If you start a
% bunch of particles at the origin on a random walk, you will find a
% correlation with diffusivity and mean-square-distance because you've
% selectively chosen the particles that diffused the furthest.
%
% The slope method is computing a slope for each pair, which presumably is
% then binned and averaged. However, it may be better to bin pair sizes
% first, then do a single slope to the whole ensemble.
%
nDrifters = size(x,2);
nReps = nDrifters*(nDrifters-1)/2;
r2 = zeros(nReps,1);
kappa_r = zeros(nReps,1);
velocity_correlation = zeros(nReps,1);
velocities = zeros(length(t)-1,nReps);
dt = t(2)-t(1);

stride = 1;
iteration = 1;
for iDrifter=1:stride:nDrifters
    for jDrifter = (iDrifter+1):stride:nDrifters        
        q = x(:,iDrifter) - x(:,jDrifter);
        r = y(:,iDrifter) - y(:,jDrifter);
        
        % mean-squared-separation distance over the requested time-interval
        r2(iteration) = mean(q.^2 + r.^2,1);
        
        % Now remove the initial conditions.
        q = q-q(1);
        r = r-r(1);
        
        % squared-separation distance as a function of time.
        D2 = q.^2 + r.^2;
        
        if strcmp(method,'endpoint')
            kappa_r(iteration) = (D2(end) - D2(1))/(4*(t(end) - t(1)));
        elseif strcmp(method,'slope')
            [p,~,mu]=polyfit(t,D2,1);
            kappa_r(iteration) = (p(1)/mu(2))/4;
        elseif strcmp(method,'powspec')
            velocities(:,iteration) = diff(q+sqrt(-1)*r)/dt;
        elseif strcmp(method,'y0v')
            % This is not a complete measure of diffusivity---this is a
            % term that we hope goes to zero. Note that if you take the
            % mean of the relative velocity (in time), that will be
            % positive, if the particles are separating on average.
            kappa_r(iteration) = q(1)*(q(end)-q(end-1))/dt + r(1)*(r(end)-r(end-1))/dt;
        elseif strcmp(method,'yv')
            % This is identical to the 'endpoint' method.
            u = diff(q)/dt;
            v = diff(r)/dt;
            qc = q(2:end) - diff(q)/2;
            rc = r(2:end) - diff(r)/2;
            kappa_r(iteration) = (mean(u.*qc) + mean(v.*rc))/2;
        end
        
        
        
%         [~,~,q,r] = CenterOfMass( x(:,[iDrifter jDrifter]), y(:,[iDrifter jDrifter]) );
        dt = t(2)-t(1);
        u = diff(x(:,[iDrifter jDrifter]))/dt;
        v = diff(y(:,[iDrifter jDrifter]))/dt;
        u = u - mean(u,1);
        v = v - mean(v,1);
        f = @(u1,u2) mean(u1.*u2)/(std(u1,1)*std(u2,1));
        velocity_correlation(iteration) = (f(u(:,1),u(:,2)) + f(v(:,1),v(:,2)))/2;
        
        iteration = iteration+1;
    end
end

if strcmp(method,'powspec')
    averaging_bandwidth = 1;
    taper_bandwidth = 0;
    kappa_r = DiffusivityFromZeroFrequency(t(2)-t(1),velocities,averaging_bandwidth,taper_bandwidth)';
end

end

