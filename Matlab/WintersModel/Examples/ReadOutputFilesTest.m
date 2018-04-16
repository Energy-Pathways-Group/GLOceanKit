% file = '/Users/jearly/Dropbox/Documents/model_raw/EarlyEtal_GM_LIN_unforced_3600000s_restart';
% file = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_unforced_36000s';
% file = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_35e-11_36000s';
file = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_NL_35e-11_36000s_restart';
% file = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/EarlyEtal_GM_LIN_unforced_3600000s_restart';

WM = WintersModel(file);
wavemodel = WM.WaveModelFromInitialConditionsFile;

nT = WM.NumberOfTimeSteps;
increments = 1:100:nT;

if exist('A_p','var') == 0
    A_p = zeros(size(increments));
    A_m = zeros(size(increments));
    B = zeros(size(increments));
    
    k = 10;
    l = 10;
    j = 1;
    
    for i = 1:length(increments)
        iTime = increments(i);
        [t,u,v,w,rho_prime] = WM.VariableFieldsAtTimeIndex(iTime,'t','u','v','w','rho_prime');
        zeta = rho_prime * wavemodel.g / (wavemodel.rho0 * wavemodel.N0);
        
        totalEnergy = mean(mean(mean( u.^2 + v.^2 + w.^2 + zeta.*zeta ) ) )/2;
        fprintf('Extracting time: %f with energy %f\n', t, totalEnergy);
        
        wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(0,u,v,rho_prime);
        
        A_p(i) = wavemodel.Amp_plus(k,l,j);
        A_m(i) = wavemodel.Amp_minus(k,l,j);
        B(i) = wavemodel.B(k,l,j);
    end
end

g = wavemodel.g;
K2 = wavemodel.K2;
h = wavemodel.h;
D = wavemodel.Lz;
f0 = wavemodel.f0;

P0_factor = sqrt(g*D*(g*h.*K2 + f0*f0)./(4*f0*f0*h)); % m^(1/2)/s
P0 = B*P0_factor(k,l,j);

E = sqrt(A_p.*conj(A_p) + A_m.*conj(A_m) + P0.*conj(P0));
figure, plot(abs(A_p)), hold on, plot(abs(A_m)), plot(abs(P0)), plot(E)