file = '/Volumes/OceanTransfer/IsotropicExperiments/QGfPlaneTurbulenceFloats_StrongForcing.mat';
output_file = '/Volumes/OceanTransfer/IsotropicExperiments/QGfPlaneTurbulenceFloatsMaternFits_Strong.mat';

load(file)

addpath('/Users/jearly/Dropbox/Documents/Projects/ForcedDissipativeQGTurbulence/MaternFit')

dt = t(2)-t(1);
numDrifters = size(cv,2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% DiffusivityFromMaternFit
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

matern_time_scale = 2*pi*dt;
frequency_scale = 1/matern_time_scale;
shouldUseZeroFrequency = 0;
shouldDisplayMaternFits = 0;
shouldAverageMaternFits = 0;
cv_dimensionless = cv; % use the same time units (as above) for the velocity
maxF = 1/(2*dt);

A_bar = zeros(size(1,numDrifters));
h_bar = zeros(size(1,numDrifters));
delta = zeros(size(1,numDrifters));
pct_fit = zeros(size(1,numDrifters));
for iDrifter=1:numDrifters
    v_mean = vmean(sqrt( cv(:,iDrifter).*conj(cv(:,iDrifter)) ), 1);
    pct_fit(iDrifter) = min(v_mean/dx/maxF,0.8);
    %     fprintf('percentage: ±%.2f\n', pct)
    fprintf('.')
    if mod(iDrifter,80)==0
        fprintf('\n');
    end
    
    [out] = MaternFit(cv_dimensionless(:,iDrifter),1/(2*pi),pct_fit(iDrifter),pct_fit(iDrifter),shouldUseZeroFrequency,shouldDisplayMaternFits,shouldAverageMaternFits);
    A_bar(iDrifter) = out.amplitude;
    h_bar(iDrifter) = out.damping;
    delta(iDrifter) = out.slope;
    
    % 	[A,d,h,omobar1]=specfit(cv_dimensionless(:,iDrifter),30,0,pi*pct,[2 0 0],psi(:,1));
    % 	A_bar(iDrifter) = A/sqrt(2*pi);
    % 	h_bar(iDrifter) = h;
    % 	delta(iDrifter) = d;
end
fprintf('\n')

% Scale to be unitful (although, fractional dimension in the case of A)
h = h_bar*frequency_scale;
A = (matern_time_scale.^(0.5-delta)) .* A_bar;

T_decorrelation = matern_time_scale./(h_bar .* beta(1/2,delta));
kappa = dt*A_bar.^2./(h_bar.^(2*delta))/4;

note = 'The matern time scale is 2*pi*dt, and has been used to make h and A unitful. I am not sure if the 2*pi should be there for A.';
save(output_file,'A_bar','h_bar', 'delta', 'A', 'h', 'T_decorrelation', 'kappa', 'matern_time_scale', 'note');

% median non-dimensional parameters
slope = median(delta);
h_typical = matern_time_scale/(median(T_decorrelation)*beta(1/2,slope));
A_typical = sqrt(4* (h_typical^(2*slope)) * median(kappa)/dt);


