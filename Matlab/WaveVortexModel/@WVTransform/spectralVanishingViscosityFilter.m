function [Qkl,Qj,kl_cutoff] = spectralVanishingViscosityFilter(self,options)
arguments
    self WVTransform {mustBeNonempty}
    options.shouldAssumeAntialiasing double {mustBeMember(options.shouldAssumeAntialiasing,[0 1])} = 1
end
% Builds the spectral vanishing viscosity operator
k_max = max(self.k);
l_max = max(self.l);
j_max = max(self.j);
if options.shouldAssumeAntialiasing == 1
    k_max = 2*k_max/3;
    l_max = 2*l_max/3;
    j_max = 2*j_max/3;
end

kl_max = min(k_max,l_max);
dkl_min = min(self.k(2)-self.k(1), self.l(2)-self.l(1));
kl_cutoff = dkl_min*(kl_max/dkl_min)^(3/4);

[K,L,J] = ndgrid(self.k,self.l,self.j);
Kh = sqrt(K.^2 + L.^2);

Qkl = exp( - ((abs(Kh)-kl_max)./(abs(Kh)-kl_cutoff)).^2 );
Qkl(abs(Kh)<kl_cutoff) = 0;
Qkl(abs(Kh)>kl_max) = 1;

if self.Nj > 2
    dj = self.j(2)-self.j(1);
    j_cutoff = dj*(j_max/dj)^(3/4);
    Qj = exp( - ((J-j_max)./(J-j_cutoff)).^2 );
    Qj(J<j_cutoff) = 0;
    Qj(J>j_cutoff) = 1;
else
    Qj = ones(size(J));
end
end