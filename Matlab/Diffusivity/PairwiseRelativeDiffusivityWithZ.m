function [r2_h, kappa_r_h, r2_z, kappa_r_z, uv_correlation, w_correlation] = PairwiseRelativeDiffusivityWithZ( t, x, y, z, method )
%PairwiseRelativeDiffusivity Compute the relative diffusivity using all
% pairwise combinations of particles.
%
% Note that the mean-square-distance, r2, can bias a sample. If you start a
% bunch of particles at the origin on a random walk, you will find a
% correlation with diffusivity and mean-square-distance because you've
% selectively chosen the particles that diffused the furthest.
nDrifters = size(x,2);
nReps = nDrifters*(nDrifters-1)/2;
r2_h = zeros(nReps,1);
kappa_r_h = zeros(nReps,1);
r2_z = zeros(nReps,1);
kappa_r_z = zeros(nReps,1);
uv_correlation = zeros(nReps,1);
w_correlation = zeros(nReps,1);
velocities_h = zeros(length(t)-1,nReps);
velocities_z = zeros(length(t)-1,nReps);
dt = t(2)-t(1);

stride = 1;
iteration = 1;
for iDrifter=1:stride:nDrifters
    for jDrifter = (iDrifter+1):stride:nDrifters        
        q = x(:,iDrifter) - x(:,jDrifter);
        r = y(:,iDrifter) - y(:,jDrifter);
        s = z(:,iDrifter) - z(:,jDrifter);
        
        % mean-squared-separation distance over the requested time-interval
        r2_h(iteration) = mean(q.^2 + r.^2,1);
        r2_z(iteration) = mean(s.^2,1);
        
        % Now remove the initial conditions.
        q = q-q(1);
        r = r-r(1);
        s = s-s(1);
        
        % squared-separation distance as a function of time.
        D2 = q.^2 + r.^2;
        D2z = s.^2;
        
        if strcmp(method,'endpoint')
            kappa_r_h(iteration) = (D2(end) - D2(1))/(4*(t(end) - t(1)));
            kappa_r_z(iteration) = (D2z(end) - D2z(1))/(2*(t(end) - t(1)));
        elseif strcmp(method,'slope')
            [p,~,mu]=polyfit(t,D2,1);
            kappa_r_h(iteration) = (p(1)/mu(2))/4;
            [p,~,mu]=polyfit(t,D2z,1);
            kappa_r_z(iteration) = (p(1)/mu(2))/2;
        elseif strcmp(method,'powspec')
            velocities_h(:,iteration) = diff(q+sqrt(-1)*r)/dt;
            velocities_z(:,iteration) = diff(s)/dt;
        elseif strcmp(method,'y0v')
            % This is not a complete measure of diffusivity---this is a
            % term that we hope goes to zero. Note that if you take the
            % mean of the relative velocity (in time), that will be
            % positive, if the particles are separating on average.
            kappa_r_h(iteration) = q(1)*(q(end)-q(end-1))/dt + r(1)*(r(end)-r(end-1))/dt;
            kappa_r_z(iteration) = s(1)*(s(end)-s(end-1))/dt;
        elseif strcmp(method,'yv')
            % This is identical to the 'endpoint' method.
            u = diff(q)/dt;
            v = diff(r)/dt;
            w = diff(s)/dt;
            qc = q(2:end) - diff(q)/2;
            rc = r(2:end) - diff(r)/2;
            sc = s(2:end) - diff(s)/2;
            kappa_r_h(iteration) = (mean(u.*qc) + mean(v.*rc))/2;
            kappa_r_z(iteration) = mean(w.*sc)/2;
        end
        
        
        
%         [~,~,q,r] = CenterOfMass( x(:,[iDrifter jDrifter]), y(:,[iDrifter jDrifter]) );
        dt = t(2)-t(1);
        u = diff(x(:,[iDrifter jDrifter]))/dt;
        v = diff(y(:,[iDrifter jDrifter]))/dt;
        w = diff(z(:,[iDrifter jDrifter]))/dt;
        u = u - mean(u,1);
        v = v - mean(v,1);
        w = w - mean(v,1);
        f = @(u1,u2) mean(u1.*u2)/(std(u1,1)*std(u2,1));
        uv_correlation(iteration) = (f(u(:,1),u(:,2)) + f(v(:,1),v(:,2)))/2;
        w_correlation(iteration) = f(w(:,1),w(:,2));
        
        iteration = iteration+1;
    end
end

if strcmp(method,'powspec')
    averaging_bandwidth = 1;
    taper_bandwidth = 0;
    kappa_r_h = DiffusivityFromZeroFrequency(t(2)-t(1),velocities_h,averaging_bandwidth,taper_bandwidth)';
    kappa_r_z = DiffusivityFromZeroFrequency(t(2)-t(1),velocities_z,averaging_bandwidth,taper_bandwidth)';
end

end

