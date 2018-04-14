file = '/Users/jearly/Dropbox/Documents/model_raw/EarlyEtal_GM_LIN_unforced_3600000s_restart';

WM = WintersModel(file);
wavemodel = WM.WaveModelFromInitialConditionsFile;

nT = WM.NumberOfTimeSteps;

iTime = 1;
[t,u,v,rho_prime] = WM.VariableFieldsAtTimeIndex(iTime,'t','u','v','rho_prime');
wavemodel.InitializeWithHorizontalVelocityAndIsopycnalDisplacementFields(t,u,v,rho_prime);

t
A1_p = wavemodel.Amp_plus;
A1_m = wavemodel.Amp_minus;
B1 = wavemodel.B;

iTime = 4;
[t,u,v,rho_prime, rho] = WM.VariableFieldsAtTimeIndex(iTime,'t','u','v','rho_prime', 'rho');
wavemodel.InitializeWithHorizontalVelocityAndIsopycnalDisplacementFields(t,u,v,rho_prime);

t
A2_p = wavemodel.Amp_plus;
A2_m = wavemodel.Amp_minus;
B2 = wavemodel.B;