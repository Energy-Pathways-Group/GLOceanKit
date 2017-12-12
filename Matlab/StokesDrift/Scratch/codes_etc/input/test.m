%program to test a nonlinear density stratification
clf
Lz=1000.
delta=0.2*Lz %governs the steepness of the transition layer
d=0.% vertical offset from z0
z0=2*Lz/3
zr=z0-d
nz=128
s1bar=zeros(nz,1);
z=zeros(nz,1);
dz=Lz/nz
for k=1:nz
z(k)=k*dz
s1bar(k)=0.5*(1.-tanh(2*(z(k)-zr)/delta))
%s1bar(k)=tanh(z(k))
end
plot(z,s1bar)
