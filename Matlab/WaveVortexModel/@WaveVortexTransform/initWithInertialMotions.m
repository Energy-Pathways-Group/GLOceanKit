function initWithInertialMotions(self,u,v)
self.Ap = zeros(size(self.Ap));
self.Am = zeros(size(self.Am));
self.A0 = zeros(size(self.A0));
self.setInertialMotions(u,v);
end