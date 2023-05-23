function du = diffZG(self,u,n)
arguments
    self        WVTransform
    u (:,:,:)   double
    n (1,1)     double = 1
end

if n ~= 1
    error('Not yet implemented')
end

u = permute(u,[3 1 2]); % keep adjacent in memory
u = reshape(u,self.Nz,[]);

D = self.PFinv*(squeeze(self.P./(self.Q .* self.h)).*self.QG);
du =  D*u;

du = reshape(du,self.Nz,self.Nx,self.Ny);
du = permute(du,[2 3 1]);

end