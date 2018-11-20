file = '/Users/jearly/Documents/ProjectRepositories/single-wave-exponential-stratification/WintersModelRuns/output_180506';

WM = WintersModel(file);

if WM.NumberOf3DOutputFiles > 1
    [t,u,v,w,x,z] = WM.VariableFieldsFrom3DOutputFileAtIndex(1,'t','u','v','w','x','z');
elseif WM.NumberOf2DOutputFiles > 1
    [t,u,v,w,x,z] = WM.VariableFieldsFrom2DOutputFileAtIndex(90,'t','u','v','w','x','z');
end

U = max(max(max(sqrt(u.*u + v.*v))));
W = max(max(max(abs(w))));

dx = x(2)-x(1);
dz = z(2)-z(1);

T_diss_x_theory = (dx/pi)^2/(U*dx);
T_diss_z_theory = (dz/pi)^2/(W*dz);

T_diss_x = WM.T_diss_x;
T_diss_z = WM.T_diss_z;

cfl = 0.25; % U*dt/dx < cfl
dT_x = cfl*dx/U;
dT_z = cfl*dz/W;

return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Compute the damping scales
%

% Kraig Winters uses an e-fold time to set \nu_x in the hypervicous
% operator. We start by computing \nu_x and \nu_z.

nu_x = (-1)^(WM.p_x+1)*power(dx/pi,2*WM.p_x) / WM.T_diss_x;
nu_z = (-1)^(WM.p_z+1)*power(dz/pi,2*WM.p_z) / WM.T_diss_z;

threshold = 0.6; % exp(-1) should get you all the way to the Nyquist.
tau = max(WM.T_diss_x,WM.T_diss_z);
m = (pi/(length(wavemodel.j)*dz))*wavemodel.j;
lambda_x = nu_x*(sqrt(-1)*wavemodel.k).^(2*WM.p_x);
lambda_z = nu_z*(sqrt(-1)*m).^(2*WM.p_z);
maxKindex = find(exp(lambda_x*tau)<threshold,1,'first');
maxJindex = find(exp(lambda_z*tau)<threshold,1,'first');

maxK = wavemodel.k(maxKindex);
maxMode = maxJindex;

nK = length(wavemodel.k)/2 + 1;
k_diss = abs(wavemodel.k(1:nK));
j_diss = 0:max(wavemodel.j); % start at 0, to help with contour drawing
[K,J] = ndgrid(k_diss,j_diss);
M = (2*pi/(length(wavemodel.j)*dz))*J/2;

lambda_x = nu_x*(sqrt(-1)*K).^(2*WM.p_x);
lambda_z = nu_z*(sqrt(-1)*M).^(2*WM.p_z);
% tau = WM.VariableFieldsFrom3DOutputFileAtIndex(WM.NumberOf3DOutputFiles,'t');
tau = max(WM.T_diss_x,WM.T_diss_z);
R = exp((lambda_x+lambda_z)*tau);

% The highest wavenumber should e-fold in time tau, so let's contour the
% area that retains 90% of its value
C = contourc(j_diss,k_diss,R,0.60*[1 1]);
n = C(2,1);
j_damp = C(1,1+1:n);
k_damp = C(2,1+1:n);