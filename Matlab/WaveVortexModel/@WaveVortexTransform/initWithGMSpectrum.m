function initWithGMSpectrum(self, GMAmplitude, varargin)
% initialize with a Garrett-Munk spectrum
%
% This only initializes the wave components, A0 is left untouched.
%
% - Topic: Initial conditions — Waves
% - Declaration: initWithGMSpectrum( GMAmplitude, varargin)
if mod(length(varargin),2) ~= 0
    error('Arguments must be given as name/value pairs.');
end

% Set defaults
j_star = 3;

% Now override the defaults with user settings
for iArg = 1:2:length(varargin)
    if strcmp(varargin{iArg}, 'j_star')
        j_star = varargin{iArg+1};
        varargin(iArg+1) = [];
        varargin(iArg) = [];
        break;
    end
end

% GM Parameters
L_gm = 1.3e3; % thermocline exponential scale, meters
invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_gm = 6.3e-5; % non-dimensional energy parameter
E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm*GMAmplitude;
%             E = E*(self.Lz/L_gm); % This correction fixes the amplitude so that the HKE variance at a given depth matches (instead of depth integrated energy)

% Compute the proper vertical function normalization
H = (j_star+(1:1024)).^(-5/2);
H_norm = 1/sum(H);

% Do the same for the frequency function.
B_norm = 1/atan(sqrt(self.Nmax*self.Nmax/(self.f0*self.f0)-1));

% This function tells you how much energy you need between two
% frequencies for a given vertical mode.
GM2D_int = @(omega0,omega1,j) E*H_norm*B_norm*((j+j_star).^(-5/2))*(atan(self.f0/sqrt(omega0*omega0-self.f0*self.f0)) - atan(self.f0/sqrt(omega1*omega1-self.f0*self.f0)));

% Do a quick check to see how much energy is lost due to
% limited vertical resolution.
maxMode = self.Nj;
for iArg = 1:2:length(varargin)
    if strcmp(varargin{iArg}, 'maxMode')
        maxMode = varargin{iArg+1};
    end
end

totalEnergy = 0;
for mode=1:maxMode
    totalEnergy = totalEnergy + GM2D_int(self.f0,self.Nmax,mode);
end
fprintf('You will miss %.2f%% of the energy due to limited vertical modes.\n',100-100*totalEnergy/E);

[GM3Dint,GM3Dext] = self.initWithSpectralFunction(GM2D_int,varargin{:});

fprintf('After distributing energy across frequency and mode, you still have %.2f%% of reference GM energy.\n',100*(sum(sum(sum(GM3Dint))) + sum(GM3Dext))/E);
fprintf('Due to restricted domain size, the j=1,k=l=0 mode contains %.2f%% the total energy.\n',100*GM3Dint(1,1,1)/(sum(sum(sum(GM3Dint))) + sum(GM3Dext)) );

GM_sum_int = sum(sum(sum(GM3Dint)))/E;
GM_sum_ext = sum(GM3Dext)/E;
C = self.Apm_TE_factor;
GM_random_sum_int = sum( C(:).*(self.Ap(:).*conj(self.Ap(:)) + self.Am(:).*conj(self.Am(:))) )/E;
GM_random_sum_ext = sum((self.offgridModes.U_cos_ext.*self.offgridModes.U_cos_ext + self.offgridModes.V_cos_ext.*self.offgridModes.V_cos_ext).*self.offgridModes.h_ext/2)/E;
fprintf('The (gridded, external) wave field sums to (%.2f%%, %.2f%%) GM given the scales, and the randomized field sums to (%.2f%%, %.2f%%) GM\n', 100*GM_sum_int, 100*GM_sum_ext, 100*GM_random_sum_int,100*GM_random_sum_ext);
end