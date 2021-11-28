function [Qk,Ql,Qj] = SpectralVanishingViscosityFilter(self)
% Builds the spectral vanishing viscosity operator
k_max = max(self.k);
l_max = max(self.l);
j_max = max(self.j);
if self.shouldAntiAlias == 1
    k_max = 2*k_max/3;
    l_max = 2*l_max/3;
    j_max = 2*j_max/3;
end

[K,L,J] = ndgrid(self.k,self.l,self.j);

dk = self.k(2)-self.k(1);
dl = self.l(2)-self.l(1);
dj = self.j(2)-self.j(1);
k_cutoff = dk*(k_max/dk)^(3/4);
l_cutoff = dl*(l_max/dl)^(3/4);
j_cutoff = dj*(j_max/dj)^(3/4);

Qk = exp( - ((abs(K)-k_max)./(abs(K)-k_cutoff)).^2 );
Qk(abs(K)<k_cutoff) = 0;
Qk(abs(K)>k_max) = 1;

Ql = exp( - ((abs(L)-l_max)./(abs(L)-l_cutoff)).^2 );
Ql(abs(L)<l_cutoff) = 0;
Ql(abs(L)>l_max) = 1;

Qj = exp( - ((J-j_max)./(J-j_cutoff)).^2 );
Qj(J<j_cutoff) = 0;
Qj(J>j_cutoff) = 1;
end