function du = diffZF(self,u,options)
arguments
    self        WVTransform
    u (:,:,:)   double
    options.n (1,1)     double = 1
end
n = options.n;
u = permute(u,[3 1 2]); % keep adjacent in memory
u = reshape(u,self.Nz,[]);

if n == 1
    DzF = -(self.N2/self.g) .* self.QG0inv*(squeeze(self.Q0 ./ self.P0).*self.PF0);
    D = DzF;
elseif n == 2
    DzF = -(self.N2/self.g) .* self.QG0inv*(squeeze(self.Q0 ./ self.P0).*self.PF0);
    DzG = self.PF0inv*(squeeze(self.P0./(self.Q0 .* self.h_0)).*self.QG0);
    D = DzG*DzF;
elseif n == 3
    DzF = -(self.N2/self.g) .* self.QG0inv*(squeeze(self.Q0 ./ self.P0).*self.PF0);
    DzzG = -(self.N2/self.g) .* self.QG0inv*(squeeze(1./self.h_0).*self.QG0);
    D = DzzG*DzF;
elseif n == 4
    DzF = -(self.N2/self.g) .* self.QG0inv*(squeeze(self.Q0 ./ self.P0).*self.PF0);
    DzzG = -(self.N2/self.g) .* self.QG0inv*(squeeze(1./self.h_0).*self.QG0);
    DzG = self.PF0inv*(squeeze(self.P0./(self.Q0 .* self.h_0)).*self.QG0);
    D = DzG * DzzG * DzF;
else
    error('Not yet implemented')
end
du =  D*u;

du = reshape(du,self.Nz,self.Nx,self.Ny);
du = permute(du,[2 3 1]);

end