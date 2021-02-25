function F = NonlinearFluxWithFloatsAndDriftersAtTimeArray(self,t,Y0,z_d)
[Fp,Fm,F0,u,v,w,ud,vd] = self.NonlinearFluxWithFloatsAndDriftersAtTime(t,Y0{1},Y0{2},Y0{3},Y0{4},Y0{5},Y0{6},Y0{7},Y0{8},z_d);
F = {Fp,Fm,F0,u,v,w,ud,vd};
end