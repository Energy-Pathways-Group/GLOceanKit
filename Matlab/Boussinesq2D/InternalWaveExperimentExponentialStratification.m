%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Adjustable parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
wavelength = 10e3; % wavelength in meters
j = 1; % vertical mode number
epsilon = 0.05; % nonlinearity parameter
maxOscillations = 5; % Total number of oscillations, in periods
stratification = 'exponential';
z0 = [-10; -250; -625]; % initial particle positions

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Derived parameters & model setup
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[rho, N2, zIn] = InternalModes.StratificationProfileWithName(stratification);
latitude = 0.0;
k = 2*pi/wavelength;
wavenumber = 4;
Nx = 128;
Nz = 513; % Must include end point to advect at the surface, so use 2^N + 1
Lx = wavenumber*wavelength;
Lz = max(zIn)-min(zIn);
z = linspace(min(zIn),max(zIn),Nz)';



model = Boussinesq2D([Lx, Lz], [Nx, Nz], z, rho);
model.nonlinear = 1;

% wwDecomp = WaveWaveDecomposition([Lx, Lz], [Nx, Nz], z, rho);
model.internalModes.normalization = Normalization.uMax;
[~,G,h,omega] = model.internalModes.ModesAtWavenumber(k);
cp = omega(j)/k;
U = epsilon*(omega(j)/k);
m = j*pi/Lz;
[omega,h] = model.InitializeWithNonlinearPlaneWave(k,j,U);
%[omega,h] = model.InitializeWithPlaneWave(k,j,U);
W = (k/m)*U;

nParticles = 5;
model.setParticlePositions(Lx*ones(1,nParticles)/2,(Lz/6)*(1:nParticles) - Lz)

% model.StepForwardToTime(5*model.dt);
% 
[A_plus,A_minus] = wwDecomp.decomposeWithWaveSpeed(model.psi,model.b,cp,model.t);



% figure
% pcolor(model.X,model.Z,model.b), shading interp, hold on
% quiver(model.X,model.Z,model.u,model.w)
% 
% model.internalModes.normalization = Normalization.kConstant;
% model.internalModes.nModes = 130;
% [~,G,h] = model.internalModes.ModesAtWavenumber(k);
% 
% psi_bar = FourierTransformForward(model.x,model.psi,1);
% b_bar = FourierTransformForward(model.x,model.b ./ model.N2,1);
% 
% psi_bar = fftshift(psi_bar,1);
% b_bar = fftshift(b_bar,1);
% % figure
% % subplot(1,2,1)
% % pcolor(abs(psi_bar)), shading flat
% % subplot(1,2,2)
% % pcolor(abs(b_bar)), shading flat
% 
% psi_decomp0 = zeros(size(psi_bar,1),size(G,2));
% b_decomp0 = zeros(size(psi_bar,1),size(G,2));
% for i=1:size(psi_bar,1)
%    psi_decomp0(i,:) =  G\(psi_bar(i,:)).';
%    b_decomp0(i,:) =  G\(b_bar(i,:)).';
% end
% figure
% subplot(1,2,1)
% pcolor(log10(abs(psi_decomp0)).'), shading flat
% subplot(1,2,2)
% pcolor(log10(abs(b_decomp0)).'), shading flat

n_period = 2*pi/omega/model.dt;
n = 10*n_period;
maxU = zeros(n,1);
maxW = zeros(n,1);
A_p = zeros(n,size(A_plus,1),size(A_plus,2));
A_m = zeros(n,size(A_plus,1),size(A_plus,2));
xi = zeros(n,nParticles);
zeta = zeros(n,nParticles);
for i=1:n
    model.StepForwardToTime((i-1)*model.dt);
    [A_plus,A_minus] = wwDecomp.decomposeWithWaveSpeed(model.psi,model.b,cp,model.t);
    A_p(i,:,:)=A_plus;
    A_m(i,:,:)=A_minus;
    xi(i,:) = model.xi;
    zeta(i,:) = model.zeta;
end

% A_p_rms = squeeze(max(A_p.*conj(A_p),[],1));
% A_p_rms = A_p_rms(1:(size(A_p_rms,1)/2),:);
% A_m_rms = squeeze(max(A_m.*conj(A_m),[],1));
% A_m_rms = A_m_rms(1:(size(A_m_rms,1)/2),:);

A_p_rms = squeeze(mean(A_p.*conj(A_p),1));
A_p_rms = A_p_rms(1:(size(A_p_rms,1)/2),:);
A_m_rms = squeeze(mean(A_m.*conj(A_m),1));
A_m_rms = A_m_rms(1:(size(A_m_rms,1)/2),:);

A_p_var = squeeze(var(A_p.*conj(A_p),1));
A_p_var = A_p_var(1:(size(A_p_var,1)/2),:);
A_m_var = squeeze(var(A_m.*conj(A_m),1));
A_m_var = A_m_var(1:(size(A_m_var,1)/2),:);

oscillatoryVariance = sum(sum((A_p_var + A_m_var)./(A_p_rms + A_m_rms + A_p_var + A_m_var)))


oscillatory_variance = sqrt(A_p_var./A_p_rms);


[a,Ip] = sort(A_p_rms(:),'descend');
[b,Im] = sort(A_m_rms(:),'descend');
names = cell(6,1);
figure
for i=1:3
   [i1,i2] = ind2sub(size(A_p_rms),Ip(i));
   plot(abs(A_p(:,i1,i2))), hold on
   names{2*(i-1)+1} = sprintf('k=%.1f k_0, j=%d',wwDecomp.k(i1)/k,i2);
   [i1,i2] = ind2sub(size(A_m_rms),Im(i));
   plot(abs(A_m(:,i1,i2))), hold on
   names{2*(i-1)+2} = sprintf('k=%.1f k_0, j=%d',wwDecomp.k(i1)/k,i2);
end
ylog
legend(names);


return
figure
plot([maxU, maxW])

figure
plot(xi,zeta)



psi_bar = FourierTransformForward(model.x,model.psi,1);
b_bar = FourierTransformForward(model.x,model.b ./ model.N2,1);

psi_bar = fftshift(psi_bar,1);
b_bar = fftshift(b_bar,1);
% figure
% subplot(1,2,1)
% pcolor(abs(psi_bar)), shading flat
% subplot(1,2,2)
% pcolor(abs(b_bar)), shading flat

psi_decomp = zeros(size(psi_bar,1),size(G,2));
b_decomp = zeros(size(psi_bar,1),size(G,2));
for i=1:size(psi_bar,1)
   psi_decomp(i,:) =  G\(psi_bar(i,:)).';
   b_decomp(i,:) =  G\(b_bar(i,:)).';
end
figure
subplot(1,2,1)
pcolor(log10(abs(psi_decomp)).'), shading flat
subplot(1,2,2)
pcolor(log10(abs(b_decomp)).'), shading flat

% I find about 130 modes for Nz = 513.
% kappa = InternalModes.ConditionNumberAsFunctionOfModeNumber(G);
% figure, plot(kappa)