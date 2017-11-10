function [ S ] = GarrettMunkHorizontalVelocitySpectrum( omega, latitude, rho, zIn, zOut, method )
%GarrettMunkHorizontalVelocitySpectrum 
%
% Returns S where size(S) = [m n] with length(zOut)=m, length(omega) = n;
%
% Rather than simply evaluate the function at the request frequencies, we
% actually integrate the function, and then divide by delta-frequency. This
% is because the function is infinite at f, but is integrable.
%
% This function does not perform well near turning frequencies (local
% buoyancy frequency).

if length(omega)>1 && ~isrow(omega)
   omega = reshape(omega,1,[]);
end

if length(zOut)>1 && ~iscolumn(zOut)
   zOut = reshape(zOut,[],1);
end

im = InternalModesWKBSpectral(rho,zIn,zOut,latitude);
N2 = im.N2;
f0 = im.f0;

% GM Parameters
j_star = 3;
L_gm = 1.3e3; % thermocline exponential scale, meters
invT_gm = 5.2e-3; % reference buoyancy frequency, radians/seconds
E_gm = 6.3e-5; % non-dimensional energy parameter
E = L_gm*L_gm*L_gm*invT_gm*invT_gm*E_gm;

H = (j_star+(1:3000)).^(-5/2);
H_norm = 1/sum(H);

% This function tells you how much energy you need between two frequencies
N_max = max(sqrt(N2));
B_norm = 1/acos(f0/N_max); % B_norm *at each depth*.
B_int = @(omega0,omega1) B_norm*(atan(f0/sqrt(omega0*omega0-f0*f0)) - atan(f0/sqrt(omega1*omega1-f0*f0)));

% Assume omega0 & omega1 >=0 and that omega1>omega0;
B = @(omega0,omega1) (omega1<f0 | omega1 > N_max)*0 + (omega0<f0 & omega1>f0)*B_int(f0,omega1) + (omega0>=f0*ones(size(zOut)) & omega1 <= N_max).*B_int(omega0,omega1) + (omega0<N_max & omega1 > N_max).*B_int(omega0,N_max);

C = @(omega) (abs(omega)<f0 | abs(omega) > N_max)*0 + (abs(omega) >= f0 & abs(omega) <= N_max)*( (1+f0/omega)*(1+f0/omega) );

dOmegaVector = diff(omega);
if any(dOmegaVector<0)
    error('omega must be strictly monotonically increasing.')
end

dOmega = unique(dOmegaVector);

if max(abs(diff(dOmega))) > 1e-7
    error('omega must be an evenly spaced grid');
end

dOmega = min( [f0/2,dOmega]);

S = zeros(length(zOut),length(omega));
for i=1:length(omega)
    Bomega = B( abs( omega(i) ) - dOmega/2, abs( omega(i) ) + dOmega/2 )/dOmega;
        
    S(:,i) = E* ( Bomega .* C(omega(i)) );
    
%     if (abs(omega(i)) > f0)
%         [F, ~, h] = im.ModesAtFrequency(omega(i));
%         j_max = find(h>0,1,'last');
%         
%         if ~isempty(j_max)
%             j_max = j_max-5; % Last mode may be bad.
%             H = H_norm*(j_star + (1:j_max)').^(-5/2);
%             Phi = sum( (F(:,1:j_max).^2) * H, 2);
%             
%             S(:,i) = S(:,i).*Phi;
%         end
%     end
    
    
end

if strcmp(method,'exact')
    nEVP = 128;
    nEVPMax = 512;
    min_j = 64; % minimum number of good modes we require
    
    im = InternalModes(rho,zIn,zOut,latitude, 'nEVP', nEVP);
    im.normalization = 'const_G_norm';
    
    [sortedOmegas, indices] = sort(abs(omega));
    for i = 1:length(sortedOmegas)
        if (sortedOmegas(i) > f0)
            
            [F, ~, h] = im.ModesAtFrequency(sortedOmegas(i));
            j_max = ceil(find(h>0,1,'last')/2);
            
            while( (isempty(j_max) || j_max < min_j) && nEVP < nEVPMax )
                nEVP = nEVP + 128;
                im = InternalModes(rho,zIn,zOut,latitude, 'nEVP', nEVP);
                im.normalization = 'const_G_norm';
                [F, ~, h] = im.ModesAtFrequency(sortedOmegas(i));
                j_max = ceil(find(h>0,1,'last')/2);
            end
            
            H = H_norm*(j_star + (1:j_max)').^(-5/2);
            Phi = sum( (F(:,1:j_max).^2) * (1./h(1:j_max) .* H), 2);
            S(:,indices(i)) = S(:,indices(i)).*Phi;
        end
    end
elseif strcmp(method,'wkb')
    im = InternalModes(rho,zIn,zOut,latitude, 'method', 'wkb');
    im.normalization = 'const_G_norm';
    
    [sortedOmegas, indices] = sort(abs(omega));
    for i = 1:length(sortedOmegas)
        if (sortedOmegas(i) > f0)
            
            [F, ~, h] = im.ModesAtFrequency(sortedOmegas(i));
            j_max = ceil(find(h>0,1,'last')/2);
                       
            H = H_norm*(j_star + (1:j_max)').^(-5/2);
            Phi = sum( (F(:,1:j_max).^2) * (1./h(1:j_max) .* H), 2);
            S(:,indices(i)) = S(:,indices(i)).*Phi;
        end
    end
    
elseif strcmp(method,'gm')
    S = S .* (sqrt(N2)/(L_gm*invT_gm));
    
    zerosMask = zeros(length(zOut),length(omega));
    for iDepth = 1:length(zOut)
        zerosMask(iDepth,:) = abs(omega) < sqrt(N2(iDepth));
    end
    S = S .* zerosMask;
end

% for the vertical structure function I think we should sort the
% frequencies, and try to solve with a small nEVP. If the number of valid
% modes starts to look small, then we up the resolution. We keep the
% resolution for higher frequencies, and repeat the check. This is good
% because high resolution is always needed as the frequencies get higher.


end

