wvt = WVTransformHydrostatic([800e3, 800e3, 1300], [64 64 30], N2=@(z) (5.2e-3)*(5.2e-3)*ones(size(z)));

%%
L_gm = 1.3e3; % thermocline exponential scale, meters
invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_gm = 6.3e-5; % non-dimensional energy parameter
E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm;

j_star = 3;
Nmax = invT_gm;

% Compute the proper vertical function normalization
H = (j_star+(1:1024)).^(-5/2);
H_norm = 1/sum(H);

% Do the same for the frequency function.
B_norm = 1/atan(sqrt(Nmax*Nmax/(wvt.f*wvt.f)-1));

H = @(j) H_norm*((j+j_star).^(-5/2));
B = @(omega) B_norm*wvt.f./(omega.*sqrt(omega.*omega-wvt.f*wvt.f));
GM = @(omega,j) E*H(j) .* B(omega);

totalEnergy = 0;
for j=1:20
    totalEnergy = totalEnergy + integral( @(omega) GM(omega,j),wvt.f,Nmax);
end
totalEnergy

wvt.initWithFrequencySpectrum(ApmSpectrum=GM);

% GM2D_int = @(omega0,omega1,j) E*B_norm*((j+j_star).^(-5/2))*(atan(self.f/sqrt(omega0*omega0-self.f*self.f)) - atan(self.f/sqrt(omega1*omega1-self.f*self.f)));