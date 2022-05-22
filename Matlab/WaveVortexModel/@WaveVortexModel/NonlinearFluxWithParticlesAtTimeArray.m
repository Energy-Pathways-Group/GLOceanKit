function F = NonlinearFluxWithParticlesAtTimeArray(self,t,Y0)
[Fp,Fm,F0,u,v,w] = self.NonlinearFluxWithParticlesAtTime(t,Y0{1},Y0{2},Y0{3},Y0{4},Y0{5},Y0{6});
F = {Fp,Fm,F0,u,v,w};
end