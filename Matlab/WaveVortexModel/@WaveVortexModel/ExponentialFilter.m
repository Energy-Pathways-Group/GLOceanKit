function [Qk,Ql,Qj] = ExponentialFilter(self)
% Builds an exponential filter
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

alpha = 15*log(10); % drop to machine precision (15 sigfigs double precision)

Q = @(x,x_max,x_cutoff) exp( - alpha*((abs(x)-x_cutoff)./(x_max-x_cutoff)).^2 );
Qk = Q(K,k_max,k_cutoff);
Qk(abs(K)< k_cutoff) = 1;
Qk(abs(K)> k_max) = 0;

Ql = Q(L,l_max,l_cutoff);
Ql(abs(L)< l_cutoff) = 1;
Ql(abs(L)> l_max) = 0;

Qj = Q(J,j_max,j_cutoff);
Qj(abs(J)< j_cutoff) = 1;
Qj(abs(J)> j_max) = 0;

end