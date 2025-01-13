function flag = isequal(self,other)
arguments
    self WVTransform
    other WVTransform
end
flag = isequal(self.shouldAntialias, other.shouldAntialias);
flag = flag & isequal(class(self),class(other));
flag = flag & isequal(self.x, other.x);
flag = flag & isequal(self.y, other.y);
flag = flag & isequal(self.z,other.z);
flag = flag & isequal(self.j, other.j);
flag = flag & isequal(self.k, other.k);
flag = flag & isequal(self.l, other.l);
flag = flag & isequal(self.t, other.t);
flag = flag & isequal(self.t0, other.t0);
flag = flag & isequal(self.rho0, other.rho0);
flag = flag & isequal(self.conjugateDimension, other.conjugateDimension);
end