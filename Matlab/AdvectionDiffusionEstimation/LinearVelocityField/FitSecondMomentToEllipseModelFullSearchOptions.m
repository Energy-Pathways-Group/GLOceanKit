function [parameters,error] = FitSecondMomentToEllipseModelFullSearchOptions( Mxx, Myy, Mxy, t, model, varargin)
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
searchAlgorithm = 'fminsearch';
expectedKappaRange = [0.1 5];
expectedSigmaRange = [1e-6 5e-5];
expectedZetaRange = [1e-7 5e-5];
divergenceMeasure = 'Bhattacharyya';
for k = 1:2:length(varargin)
    if strcmp(varargin{k}, 'searchAlgorithm')
        searchAlgorithm =  varargin{k+1};
    end
    if strcmp(varargin{k}, 'expectedKappaRange')
        expectedKappaRange = varargin{k+1};
    end
    if strcmp(varargin{k}, 'expectedSigmaRange')
        expectedSigmaRange =  varargin{k+1};
    end
    if strcmp(varargin{k}, 'expectedZetaRange')
        if any(varargin{k+1}<0)
           error('expectedZetaRange should only include positive values. The search will include negative values.')
        end
        expectedZetaRange =  varargin{k+1};
    end
    if strcmp(varargin{k}, 'divergenceMeasure')
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

% kappaScale = 0.25/20;
sScale = sigmaScale; %2.5e-6/20;
% thetaScale = (pi/180);

% How do we measure the error?
divergence = @(M_model) MomentTensorModelError(t,Mxx,Myy,Mxy,M_model(:,1),M_model(:,2),M_model(:,3),divergenceMeasure); %, 'Ellipse-Overlap' 'Bhattacharyya'

% By setting minError to Inf, these *will* be replaced.
if strcmp(model,'diffusive')
    minError = Inf; minKappa = initialKappa; minSigma = 0; minTheta = 0; minZeta = 0;
elseif strcmp(model,'strain-diffusive')
    minError = Inf; minKappa = initialKappa; minSigma = initialSigma; minTheta = initialTheta; minZeta = 0;
elseif strcmp(model,'vorticity-strain*-diffusive') % strain dominanted => strain rate must be greater than zeta
    minError = Inf; minKappa = initialKappa; minSigma = initialSigma; minTheta = initialTheta; minZeta = 0.5*sScale;
elseif strcmp(model,'vorticity*-strain-diffusive') % strain dominanted => strain rate must be greater than zeta
    minError = Inf; minKappa = initialKappa; minSigma = initialSigma; minTheta = initialTheta; minZeta = 2*sScale;
end

% Good initial values
kappa0 = initialKappa; sigma0 = initialSigma; theta0 = initialTheta; zeta0 = minZeta;

% optimization options

optimizationOptions = optimset('TolFun',1e-5,'TolX',1e-4); % , 'MaxFunEvals', 5000
% optimizationOptions = optimset('TolFun',1e-3,'TolX',1e-2,'Display','iter');


if strcmp(searchAlgorithm,'gridded')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Uses a gridded search algorithm, then 'falls through' to just do an
    % fminsearch starting from the best minimum.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    N_theta=16;
    thetaGrid = linspace(0,pi-pi/N_theta,N_theta)';
    N_kappa = 10;
    kappaGrid = 10.^linspace(log10(0.1),log10(10),N_kappa)';
    N_s = 16;
    sGrid = 10.^linspace(log10(1e-6),log10(5e-5),N_s)';
    
    if strcmp(model,'diffusive')
        for iKappa = 1:N_kappa
           error = divergence(MomentTensorEvolutionInStrainVorticityField( Mxx(1), Myy(1), Mxy(1), t, 0, 0, 0, kappaGrid(iKappa) ));
           if error < minError
               if error < minError
                   minError = error;
                   minKappa = kappaGrid(iKappa);
               end
           end
        end
    elseif strcmp(model,'strain-diffusive')
        for iKappa = 1:N_kappa
            for iS = 1:N_s
                for iTheta=1:N_theta
                    error = divergence(MomentTensorEvolutionInStrainVorticityField( Mxx(1), Myy(1), Mxy(1), t, 0, sGrid(iS), thetaGrid(iTheta), kappaGrid(iKappa) ));
                    if error < minError
                        minError = error;
                        minKappa = kappaGrid(iKappa); minSigma = sGrid(iS); minTheta = thetaGrid(iTheta);
                    end
                end
            end
        end
    elseif strcmp(model,'vorticity-strain*-diffusive')
        for iKappa = 1:N_kappa
            for iS = 1:N_s
                for iTheta=1:N_theta
                    for iZeta=(-(iS-1)):(iS-1)
                        if iZeta == 0
                            continue;
                        end
                        error = divergence(MomentTensorEvolutionInStrainVorticityField( Mxx(1), Myy(1), Mxy(1), t, sign(iZeta)*sGrid(abs(iZeta)), sGrid(iS), thetaGrid(iTheta), kappaGrid(iKappa) ));
                        if error < minError
                            minError = error;
                            minKappa = kappaGrid(iKappa); minSigma = sGrid(iS); minTheta = thetaGrid(iTheta); minZeta = sign(iZeta)*sGrid(abs(iZeta));
                        end
                    end
                end
            end
        end
    end        
    
    searchAlgorithm = 'fminsearch';
elseif strcmp(searchAlgorithm,'angle-gridded')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Uses gridded angle values to make sure we search through angle space
    % correctly. The algorithm then 'falls through' to just do an
    % fminsearch starting from the best minimum.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    if strcmp(model,'diffusive')
        % Nothing to do here---there is no angle. So fall through to
        % fminsearch.
    else
        N=8;
        thetaFixed = linspace(0,pi-pi/N,N)';
        for iTheta=1:N
            
            if strcmp(model,'strain-diffusive')
                model_error_fixed_angle = @(a) divergence(MomentTensorEvolutionInStrainVorticityField( Mxx(1), Myy(1), Mxy(1), t, 0, sScale*exp(a(1)), thetaFixed(iTheta), kappaScale*exp(a(2)) ));
                parameters_from_coeffs_fixed_angle = @(kappa,sigma,theta,zeta) [log(sigma/sScale), log(kappa/kappaScale)]; % [s, kappa_bar]
                coeffs_from_parameters_fixed_angle = @(a) [kappaScale*exp(a(2)); sScale*exp(a(1)); thetaFixed(iTheta); 0]; % [kappa, sigma, theta, zeta]
            elseif strcmp(model,'vorticity-strain*-diffusive')
                % kappa = kappaScale*exp(p_kappa); p_kappa = log(kappa/kappaScale);
                % s = sScale*exp(p_s); p_s = log(s/sScale)
                % sigma = s*cosh(alpha); s = sqrt( sigma^2 - zeta^2), so p_s = log(sqrt( sigma^2 - zeta^2)/sScale)
                % zeta = s*sinh(alpha);  alpha = atanh(zeta/sigma)
                % theta = thetaScale*p_theta
                model_error_fixed_angle = @(a) divergence(MomentTensorEvolutionInStrainVorticityField( Mxx(1), Myy(1), Mxy(1), t, sScale*exp(a(1))*sinh(a(3)), sScale*exp(a(1))*cosh(a(3)), thetaFixed(iTheta), kappaScale*exp(a(2)) ));
                parameters_from_coeffs_fixed_angle = @(kappa,sigma,theta,zeta) [log(sqrt( sigma^2 - zeta^2)/sScale), log(kappa/kappaScale), atanh(zeta/sigma)]; % [s, kappa_bar, alpha]
                coeffs_from_parameters_fixed_angle = @(a) [kappaScale*exp(a(2)); sScale*exp(a(1))*cosh(a(3)); thetaFixed(iTheta); sScale*exp(a(1))*sinh(a(3))]; % [kappa, sigma, theta, zeta]
            elseif strcmp(model,'vorticity*-strain-diffusive')
                % kappa = kappaScale*exp(p_kappa); p_kappa = log(kappa/kappaScale);
                % s = sScale*exp(p_s); p_s = log(s/sScale)
                % sigma = s*sinh(alpha); s = sqrt( zeta^2 - sigma^2), so p_s = log(sqrt( zeta^2 - sigma^2 )/sScale)
                % zeta = s*cosh(alpha);  alpha = atanh(sigma/zeta)
                % theta = thetaScale*p_theta
                model_error_fixed_angle = @(a) divergence(MomentTensorEvolutionInStrainVorticityField( Mxx(1), Myy(1), Mxy(1), t, sScale*exp(a(1))*cosh(a(3)), sScale*exp(a(1))*sinh(a(3)), thetaFixed(iTheta), kappaScale*exp(a(2)) ));
                parameters_from_coeffs_fixed_angle = @(kappa,sigma,theta,zeta) [log(sqrt( zeta^2 - sigma^2 )/sScale), log(kappa/kappaScale), atanh(sigma/zeta)]; % [s, kappa_bar, alpha]
                coeffs_from_parameters_fixed_angle = @(a) [kappaScale*exp(a(2)); sScale*exp(a(1))*sinh(a(3)); thetaFixed(iTheta); sScale*exp(a(1))*cosh(a(3))]; % [kappa, sigma, theta, zeta]
            end
            
            [a, error,~,~] = fminsearch( model_error_fixed_angle, parameters_from_coeffs_fixed_angle(kappa0,sigma0,theta0,zeta0), optimizationOptions );
            if error < minError
                minError = error;
                coeffs = coeffs_from_parameters_fixed_angle(a);
                minKappa = coeffs(1); minSigma = coeffs(2); minTheta = coeffs(3); minZeta = coeffs(4);
            end
        end
    end

    searchAlgorithm = 'fminsearch';
end  
    
if strcmp(searchAlgorithm,'fminsearch')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Uses fminsearch to find the optimal parameters
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if strcmp(model,'diffusive')
        model_error = @(a) divergence(MomentTensorEvolutionInStrainVorticityField( Mxx(1), Myy(1), Mxy(1), t, 0, 0, 0, kappaScale*exp(a(1)) ));
        parameters_from_coeffs = @(kappa,sigma,theta,zeta) log(kappa/kappaScale);
        coeffs_from_parameters = @(a) [kappaScale*exp(a(1)); 0; 0; 0]; % [kappa, sigma, theta, zeta]
        
        % Take the current minimum best guess to start the search
        [a, error,~,~] = fminsearch( model_error, parameters_from_coeffs(minKappa,minSigma,minTheta,minZeta), optimizationOptions );   
        
        if error < minError
            coeffs = coeffs_from_parameters(a);
            minError = error; minKappa = coeffs(1); minSigma = coeffs(2); minTheta = coeffs(3); minZeta = coeffs(4);
        end
    elseif strcmp(model,'strain-diffusive')
        model_error = @(a) divergence(MomentTensorEvolutionInStrainVorticityField( Mxx(1), Myy(1), Mxy(1), t, 0, sScale*exp(a(1)), thetaScale*a(2), kappaScale*exp(a(3)) ));
        parameters_from_coeffs = @(kappa,sigma,theta,zeta) [log(sigma/sScale), theta/thetaScale, log(kappa/kappaScale)]; % [s, theta_bar, kappa_bar]
        coeffs_from_parameters = @(a) [kappaScale*exp(a(3)); sScale*exp(a(1)); thetaScale*a(2); 0]; % [kappa, sigma, theta, zeta]
        
        % Take the current minimum best guess to start the search
        [a, error,~,~] = fminsearch( model_error, parameters_from_coeffs(minKappa,minSigma,minTheta,minZeta), optimizationOptions );
        
        if error < minError
            coeffs = coeffs_from_parameters(a);
            minError = error; minKappa = coeffs(1); minSigma = coeffs(2); minTheta = coeffs(3); minZeta = coeffs(4);
        end
    elseif strcmp(model,'vorticity-strain*-diffusive')
        % kappa = kappaScale*exp(p_kappa); p_kappa = log(kappa/kappaScale);
        % s = sScale*exp(p_s); p_s = log(s/sScale)
        % sigma = s*cosh(alpha); s = sqrt( sigma^2 - zeta^2), so p_s = log(sqrt( sigma^2 - zeta^2)/sScale)
        % zeta = s*sinh(alpha);  alpha = atanh(zeta/sigma)
        % theta = thetaScale*p_theta
%         model_error = @(a) divergence(MomentTensorEvolutionInStrainVorticityField( Mxx(1), Myy(1), Mxy(1), t, sScale*exp(a(1))*sinh(a(4)), sScale*exp(a(1))*cosh(a(4)), thetaScale*a(2), kappaScale*exp(a(3)) ));
%         parameters_from_coeffs = @(kappa,sigma,theta,zeta) [log(sqrt( sigma^2 - zeta^2)/sScale), theta/thetaScale, log(kappa/kappaScale), atanh(zeta/sigma)]; % [s, theta_bar, kappa_bar, alpha]
%         coeffs_from_parameters = @(a) [kappaScale*exp(a(3)); sScale*exp(a(1))*cosh(a(4)); thetaScale*a(2); sScale*exp(a(1))*sinh(a(4))]; % [kappa, sigma, theta, zeta]
%         
%         % Start with zeta at half the value of sigma to start the search
%         [a, error,~,~] = fminsearch( model_error, parameters_from_coeffs(minKappa,minSigma,minTheta,0.5*minSigma), optimizationOptions );
%         [b, b_error,~,~] = fminsearch( model_error, parameters_from_coeffs(minKappa,minSigma,minTheta,-0.5*minSigma), optimizationOptions );
%         if (b_error < error)
%             a = b;
%             error = b_error;
%         end
        
        model_error = @(a) vorticityStrainDominantedDiffusivity(divergence, Mxx(1), Myy(1), Mxy(1), t, sScale*exp(a(1)), thetaScale*a(2), kappaScale*exp(a(3)), expectedZetaRange );
        parameters_from_coeffs = @(kappa,sigma,theta,zeta) [log(sigma/sScale), theta/thetaScale, log(kappa/kappaScale)]; % [s, theta_bar, kappa_bar]
        coeffs_from_parameters = @(a) [kappaScale*exp(a(3)); sScale*exp(a(1)); thetaScale*a(2); 0]; % [kappa, sigma, theta, zeta]
        a = fminsearch( model_error, parameters_from_coeffs(minKappa,minSigma,minTheta,minZeta), optimizationOptions );
        
        coeffs = coeffs_from_parameters(a);
        minKappa = coeffs(1); minSigma = coeffs(2); minTheta = coeffs(3);
        zeta_from_coeff = @(p) sign(p)*exp(abs(p))*min(expectedZetaRange);
        coeff_from_zeta = @(zeta) sign(zeta)*log(abs(zeta/min(expectedZetaRange)));
        model_error = @(p_zeta) divergence(MomentTensorEvolutionInStrainVorticityField( Mxx(1), Myy(1), Mxy(1), t, zeta_from_coeff(p_zeta),minSigma, minTheta, minKappa ));
        [minP_zeta, minError] = fminbnd(model_error,coeff_from_zeta(-0.99*minSigma),coeff_from_zeta(0.99*minSigma));
        minZeta = zeta_from_coeff(minP_zeta);
    elseif strcmp(model,'vorticity*-strain-diffusive')
        % kappa = kappaScale*exp(p_kappa); p_kappa = log(kappa/kappaScale);
        % s = sScale*exp(p_s); p_s = log(s/sScale)
        % sigma = s*abs(sinh(alpha)); s = sqrt( zeta^2 - sigma^2), so p_s = log(sqrt( zeta^2 - sigma^2 )/sScale)
        % zeta = s*cosh(alpha);  alpha = atanh(sigma/zeta)
        % theta = thetaScale*p_theta
%         model_error_pos = @(a) divergence(MomentTensorEvolutionInStrainVorticityField( Mxx(1), Myy(1), Mxy(1), t, sScale*exp(a(1))*cosh(a(4)), sScale*exp(a(1))*abs(sinh(a(4))), thetaScale*a(2), kappaScale*exp(a(3)) ));
%         parameters_from_coeffs_pos = @(kappa,sigma,theta,zeta) [log(sqrt( zeta^2 - sigma^2 )/sScale), theta/thetaScale, log(kappa/kappaScale), abs(atanh(sigma/zeta))]; % [s, theta_bar, kappa_bar, alpha]
%         coeffs_from_parameters_pos = @(a) [kappaScale*exp(a(3)); sScale*exp(a(1))*abs(sinh(a(4))); thetaScale*a(2); sScale*exp(a(1))*cosh(a(4))]; % [kappa, sigma, theta, zeta]
%         
%         % Start with sigma at half the value of zeta to start the search
%         [a_pos, error_pos,~,~] = fminsearch( model_error_pos, parameters_from_coeffs_pos(minKappa,0.5*minSigma,minTheta,minSigma), optimizationOptions );
%         
%         model_error_neg = @(a) divergence(MomentTensorEvolutionInStrainVorticityField( Mxx(1), Myy(1), Mxy(1), t, -sScale*exp(a(1))*cosh(a(4)), sScale*exp(a(1))*abs(sinh(a(4))), thetaScale*a(2), kappaScale*exp(a(3)) ));
%         parameters_from_coeffs_neg = @(kappa,sigma,theta,zeta) [log(sqrt( zeta^2 - sigma^2 )/sScale), theta/thetaScale, log(kappa/kappaScale), abs(atanh(sigma/zeta))]; % [s, theta_bar, kappa_bar, alpha]
%         coeffs_from_parameters_neg = @(a) [kappaScale*exp(a(3)); sScale*exp(a(1))*abs(sinh(a(4))); thetaScale*a(2); -sScale*exp(a(1))*cosh(a(4))]; % [kappa, sigma, theta, zeta]
%         
%         [a_neg, error_neg,~,~] = fminsearch( model_error_neg, parameters_from_coeffs_neg(minKappa,0.5*minSigma,minTheta,-minSigma), optimizationOptions );
%         
%         if (error_pos < error_neg)
%             a = a_pos;
%             error = error_pos;
%             coeffs_from_parameters = coeffs_from_parameters_pos;
%         else
%             a = a_neg;
%             error = error_neg;
%             coeffs_from_parameters = coeffs_from_parameters_neg;
%         end
%         
%         if error < minError
%             coeffs = coeffs_from_parameters(a);
%             minError = error; minKappa = coeffs(1); minSigma = coeffs(2); minTheta = coeffs(3); minZeta = coeffs(4);
%         end

        model_error = @(a) vorticityStrainDiffusivity(divergence, Mxx(1), Myy(1), Mxy(1), t, sScale*exp(a(1)), thetaScale*a(2), kappaScale*exp(a(3)), expectedZetaRange );
        parameters_from_coeffs = @(kappa,sigma,theta,zeta) [log(sigma/sScale), theta/thetaScale, log(kappa/kappaScale)]; % [s, theta_bar, kappa_bar]
        coeffs_from_parameters = @(a) [kappaScale*exp(a(3)); sScale*exp(a(1)); thetaScale*a(2); 0]; % [kappa, sigma, theta, zeta]
        a = fminsearch( model_error, parameters_from_coeffs(minKappa,minSigma,minTheta,minZeta), optimizationOptions );
        
        coeffs = coeffs_from_parameters(a);
        minKappa = coeffs(1); minSigma = coeffs(2); minTheta = coeffs(3);
        zeta_from_coeff = @(p) sign(p)*exp(abs(p))*min(expectedZetaRange);
        coeff_from_zeta = @(zeta) sign(zeta)*log(abs(zeta/min(expectedZetaRange)));
        model_error = @(p_zeta) divergence(MomentTensorEvolutionInStrainVorticityField( Mxx(1), Myy(1), Mxy(1), t, zeta_from_coeff(p_zeta),minSigma, minTheta, minKappa ));
        [minP_zeta, minError] = fminbnd(model_error,coeff_from_zeta(-max(expectedZetaRange)),coeff_from_zeta(max(expectedZetaRange)));
        minZeta = zeta_from_coeff(minP_zeta);
    end
        
    error = minError;
    parameters.sigma = minSigma;
    parameters.theta = minTheta;
    parameters.kappa = minKappa;
    parameters.zeta = minZeta;
end

end

function error = vorticityStrainDominantedDiffusivity(divergence, Mxx, Myy, Mxy, t, sigma, theta, kappa, expectedZetaRange )
% Given fixed values of kappa, sigma, theta, we search across zeta using
% fminbnd to stay constrained between [-sigma, sigma].
    zeta_from_coeff = @(p) sign(p)*exp(abs(p))*min(expectedZetaRange);
    coeff_from_zeta = @(zeta) sign(zeta)*log(abs(zeta/min(expectedZetaRange)));
    % This logarithmically spans from -zeta to +zeta, never going below
    % min(zetaRange)
    model_error = @(p_zeta) divergence(MomentTensorEvolutionInStrainVorticityField( Mxx, Myy, Mxy, t, zeta_from_coeff(p_zeta) ,sigma, theta, kappa ));
    [~, error] = fminbnd(model_error,coeff_from_zeta(-0.99*sigma),coeff_from_zeta(0.99*sigma));
end

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

% 
% function error = vorticityDominatedStrainDiffusivity(divergence, Mxx, Myy, Mxy, t, zeta, sigma, theta, kappa )
% % Given fixed zeta, search across kappa, sigma, theta
%     model_error = @(zeta) divergence(MomentTensorEvolutionInStrainVorticityField( Mxx, Myy, Mxy, t, zeta,sigma, theta, kappa ));
%     
%             model_error = @(a) divergence(MomentTensorEvolutionInStrainVorticityField( Mxx, Myy, Mxy, t, zeta, sScale*exp(a(1)), thetaScale*a(2), kappaScale*exp(a(3)) ));
%         parameters_from_coeffs = @(kappa,sigma,theta,zeta) [log(sigma/sScale), theta/thetaScale, log(kappa/kappaScale)]; % [s, theta_bar, kappa_bar]
%         coeffs_from_parameters = @(a) [kappaScale*exp(a(3)); sScale*exp(a(1)); thetaScale*a(2); 0]; % [kappa, sigma, theta, zeta]
%         
%         % Take the current minimum best guess to start the search
%         [a, error,~,~] = fminsearch( model_error, parameters_from_coeffs(minKappa,minSigma,minTheta,minZeta), optimizationOptions );
%     
%     [~, error] = fminbnd(model_error,-0.99*sigma,0.99*sigma);
% end

