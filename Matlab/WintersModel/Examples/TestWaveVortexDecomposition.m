runtype = 'nonlinear';
ReadOverNetwork = 0;

if ReadOverNetwork == 1
    baseURL = '/Volumes/seattle_data1/cwortham/research/nsf_iwv/model_raw/';
else
    baseURL = '/Volumes/Samsung_T5/nsf_iwv/2018_12/';
end

if strcmp(runtype,'linear')
    file = strcat(baseURL,'EarlyV2_GM_LIN_unforced_damped_restart');
elseif strcmp(runtype,'nonlinear')
    file = strcat(baseURL,'EarlyV2_GM_NL_forced_damped_restart'); 
else
    error('invalid run type.');
end


WM = WintersModel(file);
wavemodel = WM.wavemodel;

nFiles = WM.NumberOf3DOutputFiles;
[t,u,v,w,rho_prime] = WM.VariableFieldsFrom3DOutputFileAtIndex(1,'t','u','v','w','rho_prime');
wavemodel.InitializeWithHorizontalVelocityAndDensityPerturbationFields(t,u,v,rho_prime);

[u_back,v_back,w_back,rho_prime_back] = wavemodel.VariableFieldsAtTime(t,'u','v','w','rho_prime');

error2 = @(u,u_unit) abs((u-u_unit))./(max(max(max(abs(u_unit)))));

u_error = error2(u,u_back);
max(u_error(:))

v_error = error2(v,v_back);
max(v_error(:))

w_error = error2(w,w_back);
max(w_error(:))

rho_error = error2(rho_prime,rho_prime_back);
max(rho_error(:))

% Best way to find the errors.

du = u-u_back;
dv = v-v_back;
dw = w-w_back;
drho = rho_prime - rho_prime_back;

% Currently only the nyquist frequency is bad in x
du_bar = CosineTransformForward(wavemodel.z, FourierTransformForward(wavemodel.y,FourierTransformForward(wavemodel.x,du,1),2));
[~,indices] = sort(abs(du_bar(:)),'descend');
[i,j,k] = ind2sub(size(du_bar),indices(1));

% same here---bad nyquist
dv_bar = CosineTransformForward(wavemodel.z, FourierTransformForward(wavemodel.y,FourierTransformForward(wavemodel.x,dv,1),2));
[~,indices] = sort(abs(dv_bar(:)),'descend');
[i,j,k] = ind2sub(size(dv_bar),indices(1));

% again, bad nyquist
dw_bar = SineTransformForward(wavemodel.z, FourierTransformForward(wavemodel.y,FourierTransformForward(wavemodel.x,dw,1),2));
[~,indices] = sort(abs(dw_bar(:)),'descend');
[i,j,k] = ind2sub(size(dw_bar),indices(1));

% vertical mode 0..4 etc, only at k_n=0,l_n=0
% so, we have modifications to the mean density profile. Why?!?!?!?
drho_bar = SineTransformForward(wavemodel.z, FourierTransformForward(wavemodel.y,FourierTransformForward(wavemodel.x,drho,1),2));
[~,indices] = sort(abs(drho_bar(:)),'descend');
[i,j,k] = ind2sub(size(drho_bar),indices(1));
% figure, plot(abs(squeeze(drho_bar(1,1,:))))
% figure, plot(squeeze(mean(mean(rho_prime,1),2)),wavemodel.z)

[~,index] = max(abs(dv_bar(:)));
[i,j,k] = ind2sub(size(dv_bar),index);
[rel_size,indices] = sort(abs(dv_bar(:)),'descend');