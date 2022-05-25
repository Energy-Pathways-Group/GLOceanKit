function du = diffZG(self,u,n)
arguments
    self        WaveVortexTransform
    u (:,:,:)   double
    n (1,1)     double = 1
end

u = permute(u,[3 1 2]); % keep adjacent in memory
u = reshape(u,self.Nz,[]);

% sine goes to, [1,-1,-1,1] for numDerivs = [1,2,3,4]
thesign = [1,-1,-1,1];
m = reshape(pi*self.j/self.Lz,[],1);
du_bar = thesign(mod(n-1,4)+1)*(m.^n) .* (self.DST*u);

if mod(numDerivs,2) == 0
    du = self.iDST*du_bar;
else
    du = self.iDCT*du_bar;
end
du = reshape(du,self.Nz,self.Nx,self.Ny);
du = permute(du,[2 3 1]);

end