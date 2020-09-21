function [parameters,err] = FitSecondMomentToEllipseModel( Mxx, Myy, Mxy, t, model, varargin)
% The angle-gridded search algorithm appears to be the most robust, with
% some simple testing.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Start by overriding the default options.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if mod(length(varargin),2) ~= 0
    error('Arguments must be given as name/value pairs');
end
expectedKappaRange = [0.1 5];
expectedSigmaRange = [1e-6 5e-5];
expectedZetaRange = [1e-7 5e-5];
divergenceMeasure = 'Bhattacharyya';
for k = 1:2:length(varargin)
    if strcmp(varargin{k}, 'expectedKappaRange')
        expectedKappaRange = varargin{k+1};
    elseif strcmp(varargin{k}, 'expectedSigmaRange')
        expectedSigmaRange =  varargin{k+1};
    elseif strcmp(varargin{k}, 'expectedZetaRange')
        if any(varargin{k+1}<0)
           error('expectedZetaRange should only include positive values. The search will include negative values.')
        end
        expectedZetaRange =  varargin{k+1};
    elseif strcmp(varargin{k}, 'divergenceMeasure')
        divergenceMeasure =  varargin{k+1};
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Rescale parameters to take good step sizes
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% fminsearch does not allow us to specify the initial step size,
% unfortunately. This means that we need to scale the parameters so that it
% takes an appropriate step. According to the documentation the first step
% is 5% of the initial value.
%
% Both kappa and sigma search parameters (param) are defined as follows:
%   kappa = kappaScale*exp(param);
%   log(kappa/kappaScale) = param
% The initial step will therefore be:
%   dparam0 = 0.05*log(kappa0/kappaScale)
% I know I want dparam0 to be about 0.5, and I know kappa0. So,
%   kappaScale = kappa0/exp(20*dparam0)
initialDeltaLogKappa = (log(max(expectedKappaRange))-log(min(expectedKappaRange)))/11;
initialKappa = exp((log(max(expectedKappaRange))+log(min(expectedKappaRange)))/2);
kappaScale = initialKappa/exp(20*initialDeltaLogKappa);

initialDeltaLogSigma = (log(max(expectedSigmaRange))-log(min(expectedSigmaRange)))/11;
initialSigma = exp((log(max(expectedSigmaRange))+log(min(expectedSigmaRange)))/2);
sigmaScale = initialSigma/exp(20*initialDeltaLogSigma);

% Theta has a similar problem. We need a step size of around 25 degrees. We
% define theta so that,
%   theta = thetaScale*param
%   param = theta/thetaScale
% The initial step will therefore be,
%   dparam0 = 0.05*theta0/thetaScale
%   dtheta0 = 0.05*theta0
% which means that scaling won't work. And therefore we should simply rely
% on periodicity to help us out.
%   theta0 = 20*dtheta0;

initialDeltaTheta = 25*pi/180;
initialTheta = 20*initialDeltaTheta;
thetaScale = 1;

% How do we measure the error?
divergence = @(M_model) MomentTensorModelError(t,Mxx,Myy,Mxy,M_model(:,1),M_model(:,2),M_model(:,3),divergenceMeasure); %, 'Ellipse-Overlap' 'Bhattacharyya'

optimizationOptions = optimset('TolFun',1e-5,'TolX',1e-4); % , 'MaxFunEvals', 5000
% optimizationOptions = optimset('TolFun',1e-3,'TolX',1e-2,'Display','iter');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Uses fminsearch to find the optimal parameters
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmp(model,'diffusive')
    % Define the model and the scalings
    model_error = @(a) divergence(MomentTensorEvolutionInStrainVorticityField( Mxx(1), Myy(1), Mxy(1), t, 0, 0, 0, kappaScale*exp(a(1)) ));
    parameters_from_coeffs = @(kappa,sigma,theta,zeta) log(kappa/kappaScale);
    coeffs_from_parameters = @(a) [kappaScale*exp(a(1)); 0; 0; 0]; % [kappa, sigma, theta, zeta]
    
    [a, err,~,~] = fminsearch( model_error, parameters_from_coeffs(initialKappa,initialSigma,initialTheta,0), optimizationOptions );
    
    coeffs = coeffs_from_parameters(a);
    parameters.kappa = coeffs(1); parameters.sigma = coeffs(2); parameters.theta = coeffs(3); parameters.zeta = coeffs(4);
    
elseif strcmp(model,'strain-diffusive')
    % Define the model and the scalings
    model_error = @(a) divergence(MomentTensorEvolutionInStrainVorticityField( Mxx(1), Myy(1), Mxy(1), t, 0, sigmaScale*exp(a(1)), thetaScale*a(2), kappaScale*exp(a(3)) ));
    parameters_from_coeffs = @(kappa,sigma,theta,zeta) [log(sigma/sigmaScale), theta/thetaScale, log(kappa/kappaScale)]; % [s, theta_bar, kappa_bar]
    coeffs_from_parameters = @(a) [kappaScale*exp(a(3)); sigmaScale*exp(a(1)); thetaScale*a(2); 0]; % [kappa, sigma, theta, zeta]
    
    [a, err,~,~] = fminsearch( model_error, parameters_from_coeffs(initialKappa,initialSigma,initialTheta,0), optimizationOptions );
    
    coeffs = coeffs_from_parameters(a);
    parameters.kappa = coeffs(1); parameters.sigma = coeffs(2); parameters.theta = coeffs(3); parameters.zeta = coeffs(4);
    
elseif strcmp(model,'vorticity-strain-diffusive')
    % We search *exactly* as we did over the strain-diffusive model, but
    % for each value, we do an fminbnd search over valid zeta.
    model_error = @(a) vorticityStrainDiffusivity(divergence, Mxx(1), Myy(1), Mxy(1), t, sigmaScale*exp(a(1)), thetaScale*a(2), kappaScale*exp(a(3)), expectedZetaRange );
    parameters_from_coeffs = @(kappa,sigma,theta,zeta) [log(sigma/sigmaScale), theta/thetaScale, log(kappa/kappaScale)]; % [s, theta_bar, kappa_bar]
    coeffs_from_parameters = @(a) [kappaScale*exp(a(3)); sigmaScale*exp(a(1)); thetaScale*a(2); 0]; % [kappa, sigma, theta, zeta]
    a = fminsearch( model_error, parameters_from_coeffs(initialKappa,initialSigma,initialTheta,0), optimizationOptions );
    
    % While it's true we found the minimum, this technique didn't allow us
    % to get back zeta as a returned parameter.
    coeffs = coeffs_from_parameters(a);
    parameters.kappa = coeffs(1); parameters.sigma = coeffs(2); parameters.theta = coeffs(3);
    zeta_from_coeff = @(p) sign(p)*exp(abs(p))*min(expectedZetaRange);
    coeff_from_zeta = @(zeta) sign(zeta)*log(abs(zeta/min(expectedZetaRange)));
    model_error = @(p_zeta) divergence(MomentTensorEvolutionInStrainVorticityField( Mxx(1), Myy(1), Mxy(1), t, zeta_from_coeff(p_zeta),parameters.sigma, parameters.theta, parameters.kappa ));
    [minP_zeta, err] = fminbnd(model_error,coeff_from_zeta(-max(expectedZetaRange)),coeff_from_zeta(max(expectedZetaRange)));
    parameters.zeta = zeta_from_coeff(minP_zeta);
end

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Called during the 'vorticity-strain-diffusive' search above with new
% parameters for sigma, theta kappa. This function then search within a
% range of zeta.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function error = vorticityStrainDiffusivity(divergence, Mxx, Myy, Mxy, t, sigma, theta, kappa, expectedZetaRange )
% Given fixed values of kappa, sigma, theta, we search across zeta using
% fminbnd to stay constrained between [-sigma, sigma].
    zeta_from_coeff = @(p) sign(p)*exp(abs(p))*min(expectedZetaRange);
    coeff_from_zeta = @(zeta) sign(zeta)*log(abs(zeta/min(expectedZetaRange)));
    % This logarithmically spans from -zeta to +zeta, never going below
    % min(zetaRange)
    model_error = @(p_zeta) divergence(MomentTensorEvolutionInStrainVorticityField( Mxx, Myy, Mxy, t, zeta_from_coeff(p_zeta) ,sigma, theta, kappa ));
    [~, error] = fminbnd(model_error,coeff_from_zeta(-max(expectedZetaRange)),coeff_from_zeta(max(expectedZetaRange)));
end

