% kappa = kappaScale*exp(kappa);
% s = sScale*exp(s);
% sigma = s*cosh(alpha);
% zeta = s*sinh(alpha);
% theta = theta
function [error] = MomentTensorModelError( Mxx, Myy, Mxy, t, alpha, s, theta, kappa, sScale, kappaScale, divergence )

if (kappaScale ~= 0)
    kappa = kappaScale*exp(kappa);
end
if (sScale ~= 0)
    s = sScale*exp(s);
    sigma = s*cosh(alpha);
    zeta = s*sinh(alpha);
else
    sigma = 0;
    zeta = 0;
end

[Mxx_model, Myy_model, Mxy_model] = MomentTensorEvolutionInStrainVorticityField( Mxx(1), Myy(1), Mxy(1), t, zeta, sigma, theta, kappa );

error_i = zeros(size(t));
area_i_overlap = zeros(size(t));
area_i_obs = zeros(size(t));
area_i_model = zeros(size(t));

% if ( strcmp(divergence, 'Least-Squares') )
% %     error_i = ((Mxx_model-Mxx).^2)./Mxx + ((Myy_model-Myy).^2)./Myy + ((Mxy_model-Mxy).^2)./Mxy;
%     error_i = ((Mxx_model-Mxx))./Mxx + ((Myy_model-Myy))./Myy;% + ((Mxy_model-Mxy).^2)./Mxy.^2;
% else
    for iTime=1:length(t)
        Sigma_1 = [Mxx(iTime), Mxy(iTime); Mxy(iTime), Myy(iTime)];
        Sigma_2 = [Mxx_model(iTime), Mxy_model(iTime); Mxy_model(iTime), Myy_model(iTime)];
        
        if ( strcmp(divergence, 'Ellipse-Overlap') )
            [a,b,angle] = ellipseComponentsFromMatrixComponents(Mxx(iTime),Myy(iTime),Mxy(iTime));
            [aModel,bModel,angleModel] = ellipseComponentsFromMatrixComponents(Mxx_model(iTime),Myy_model(iTime),Mxy_model(iTime));
            [overlapArea,obsArea,modelArea] = ellipse_ellipse_overlap(a,b,angle,aModel,bModel,angleModel);
            error_i(iTime) = (modelArea - overlapArea) + (obsArea - overlapArea);
            area_i_overlap(iTime) = overlapArea;
            area_i_obs(iTime) = obsArea;
            area_i_model(iTime) = modelArea;
        elseif ( strcmp(divergence, 'Symmetric-KL') )
            error_i(iTime) = 0.5*trace( Sigma_1\Sigma_2 + Sigma_2\Sigma_1 - 2*eye(2) );
        elseif ( strcmp(divergence, 'Hellinger') )
            Gamma = det( (Sigma_1 + Sigma_2)./2 );
            rho = sqrt(det(Sigma_1)*det(Sigma_2)/(Gamma*Gamma));
            error_i(iTime) = sqrt( 2*(1-rho) );
        elseif ( strcmp(divergence, 'Bhattacharyya') )
            Gamma = det( (Sigma_1 + Sigma_2)./2 );
            rho = sqrt(det(Sigma_1)*det(Sigma_2)/(Gamma*Gamma));
            error_i(iTime) = -log(rho);
        elseif ( strcmp(divergence, 'Riemann') )
            lambda = eig(Sigma_1, Sigma_2);
            error_i(iTime) = sqrt(sum(log(lambda).^2));
        elseif ( strcmp(divergence, 'Area') )
            lambda = eig(Sigma_1, Sigma_2);
            error_i(iTime) = 2 - sum((min(lambda,1./lambda)));
        end
    end
% end
if ( strcmp(divergence, 'Ellipse-Overlap') )
    error = sum(error_i)/sum(area_i_obs);
%     error = error_i(end)/area_i_obs(end);
%     fprintf('(%.2g,%.2f) ',kappa,error);
else
    error = sum(error_i);
end

%error = sum( ((Mxx_model-Mxx).^2)/sigma2_xx + ((Myy_model-Myy).^2)/sigma2_yy + ((Mxy_model-Mxy).^2)/sigma2_xy );
