function [kappa] = DiffusivityFromZeroFrequency( dt, cv, bandwidth, taper_amp )
%DiffusivityFromZeroFrequency  Compute the diffusivity by band averaging around the zero frequency
%	dt is the time interval
%	cv is the complex speed (u+iv) in meters/second
%	bandwidth indicates the number of frequencies around zero to average over
%	taper_amp indicates how aggressively you want to taper

if (taper_amp == 0)
	[f, st] = powspec(dt, cv);
	
	pos_indices = find(f>=0);
	%neg_indices = 1:(min(pos_indices)-1);
	neg_indices=( find( f>=-max(f) & f<0));
	
	f = f(pos_indices);
	
	% Note that we are *not* trying to make a one-sided spectrum here
	% We average the positive and negative components, not add them.
	s = st( pos_indices,: );
	s(2:end,:) = 0.5*(s(2:end,:) + flipud(st( neg_indices,: )));
else
	[psi,lambda]=sleptap(size(cv,1),taper_amp);
	[f,spp,snn,spn]=mspec(dt, cv,psi);
	
	% Note that we are *not* trying to make a one-sided spectrum here
	% We average the positive and negative components, not add them.
	s = spp+snn;
	s(2:end) = s(2:end)/2;
end

df = f(2)-f(1);
indices = find( f <= bandwidth*df );

% The first factor of 1/2 comes from the definition
% The second factor of 1/2 comes from the two-dimensionality.
kappa = (vmean( s(indices,:), 1))/4;
