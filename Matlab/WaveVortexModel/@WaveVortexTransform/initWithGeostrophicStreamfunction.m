function initWithGeostrophicStreamfunction(self,psi)
self.Ap = zeros(size(self.Ap));
self.Am = zeros(size(self.Am));
self.A0 = zeros(size(self.A0));
self.setGeostrophicStreamfunction(psi);
end