%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% InternalWaveModelPlaneWaveUnitTest
%
% This script uses the InternalWaveModel to create, and validate, a single
% internal wave for all possible wavenumbers.
%
% Jeffrey J. Early
% jeffrey@jeffreyearly.com
%
% November 17th, 2016   Version 1.0


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Specify the problem dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Lx = 15e3;
Ly = 15e3;
Lz = 5000;

Nx = 4;
Ny = 8;
Nz = 5;

latitude = [31]; % usually we check latitude 0 as well.
isHydrostatic = 1;
N0 = 5.2e-3/2; % Choose your stratification 7.6001e-04
U = 0.01; % m/s
phi = pi*0.232; % random, just to see if its doing the right thing API = 1 will fail, because you can't set the phase using that API.
t = 2.13*86400;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
totalErrors = 0;
totalTests = 0;

for iLat = 1:length(latitude)
    fprintf('\nlatitude: %.1f\n',latitude(iLat));
    rho0 = 1025;
    if 1 == 0
        wavemodel = WVTransformConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], N0=N0,latitude=latitude(iLat), isHydrostatic=isHydrostatic,rho0=rho0);
        rho0 = wavemodel.rho0;
        g = 9.81;
    else
        fprintf('WARNING: computing the off-grid (external) modes is slow because it has to solve the EVP each time, but do not worry, the internal transforms remain fast.\n')
        isHydrostatic = 1;
        g = 9.81;
        rho = @(z) -(N0*N0*rho0/g)*z + rho0;
        wavemodel = WVTransformHydrostatic([Lx, Ly, Lz], [Nx, Ny, Nz],latitude=latitude(iLat),N2=@(z) N0*N0*ones(size(z)));
    end
    
    % pull out some model constants for easy reference
    f = wavemodel.f;
    k = wavemodel.k;
    l = wavemodel.l;
    h = wavemodel.h_pm;
    x = wavemodel.x;
    y = wavemodel.y;
    z = wavemodel.z;
    [X,Y,Z] = ndgrid(x,y,z);
    
    for k_loop=(-Nx/2 + 1):1:(Nx/2-1)
        for l_loop=(-Ny/2 + 1):1:(Ny/2-1)
            fprintf('(k0,l0)=(%d,%d) ',k_loop,l_loop);
            for j0=0:(Nz-2)
                for thesign=-1:2:1
                    
                    if j0==0 && (k_loop ~=0 || l_loop ~= 0)
                        continue;
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %
                    % Create the analytical solution
                    %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    if (k_loop < 0)
                        k0 = Nx + k_loop;
                    else
                        k0 = k_loop;
                    end
                    if (l_loop < 0)
                        l0 = Ny + l_loop;
                    else
                        l0 = l_loop;
                    end
                    
                    m = j0*pi/Lz;
                    kk = k(k0+1);
                    ll = l(l0+1);
                    alpha=atan2(ll,kk);
                    K = sqrt( kk^2 + ll^2);
                    if j0 == 0
                        omega = f;
                    else
                        if isHydrostatic ==1
                            omega = thesign*sqrt( (K*K*N0*N0 + m*m*f*f)/(m*m) );
                        else
                            omega = thesign*sqrt( (K*K*N0*N0 + m*m*f*f)/(K*K+m*m) );
                        end
                    end
                    
                    if latitude(iLat) == 0
                        f0OverOmega = 0;
                        if K == 0
                            kOverOmega = 0;
                        else
                            kOverOmega = K/omega;
                        end
                    else
                        f0OverOmega = f/omega;
                        kOverOmega = K/omega;
                    end
                    
                    theta = kk*X + ll*Y + omega*t + phi;
                    u_unit = U*(cos(alpha)*cos( theta ) + f0OverOmega*sin(alpha)*sin( theta )).*cos(m*Z);
                    v_unit = U*(sin(alpha)*cos( theta ) - f0OverOmega*cos(alpha)*sin( theta )).*cos(m*Z);
                    w_unit = (U*K/m) * sin( theta ) .* sin(m*Z);
                    if j0 == 0
                        zeta_unit = zeros(size(Z));
                    else
                        zeta_unit = -(U/m) * kOverOmega * cos( theta ) .* sin(m*Z);
                    end
                    rho_prime_unit = -(rho0/g)*N0*N0*(U*K/m/omega) * cos( theta ) .* sin(m*Z);
                    rho_unit = rho0 -(N0*N0*rho0/g)*reshape(wavemodel.z,1,1,[]) + rho_prime_unit;
                    
                    if any(any(any(isnan(zeta_unit))))
                        error('nan')
                    end
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    %
                    % Now test how the different APIs perform at creating the
                    % requested wave.
                    %
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    for API = 1:2
                        wavemodel.removeAllExternalWaves();
                        wavemodel.removeAllWaves();
                        
                        if API == 1
                            apiName = 'initWithWaveModes';
                            wavemodel.initWithWaveModes(k=k_loop,l=l_loop,j=j0,phi=phi,U=U,sign=thesign);         
                        elseif API == 2
                            apiName = 'setWaveModes';
                            wavemodel.setWaveModes(k=k_loop,l=l_loop,j=j0,phi=phi,U=U,sign=thesign);
                            [omega_p, alpha_p, k_p, l_p, j_p, phi_p, A_p,norm] = wavemodel.waveModesFromWaveCoefficients();
                        elseif API == 3
                            apiName = 'setExternalWavesWithWavenumbers';
                            if omega_p == 0
                                signchange = 1;
                            else
                                signchange = sign(omega_p);
                            end
                            if ~isempty(k_p)
                                wavemodel.setExternalWavesWithWavenumbers(signchange*k_p, signchange*l_p, j_p, signchange*phi_p,signchange*A_p,norm);
                            end
                        elseif API == 4
                            apiName = 'setExternalWavesWithFrequencies';
                            if ~isempty(k_p)
                                wavemodel.setExternalWavesWithFrequencies(omega_p, alpha_p, j_p, phi_p,A_p,norm);
                            end
                        end
                        
%                         [u,v,w,zeta,rho_prime] = wavemodel.VariableFieldsAtTime(t,'u','v','w','eta','rho_prime');

wavemodel.t = t;
[u,v,w,zeta,rho_prime] = wavemodel.variableWithName('u','v','w','eta','rho_prime');

%                         rho = wavemodel.DensityFieldAtTime(t);
                        
                        % Compute the relative error
                        max_speed = max(max(max( sqrt(u.*u + v.*v) )));
                        max_u = max( [max(max(max( u ))), 1e-15] );
                        u_error = max(max(max(abs(u-u_unit)/max_speed)));
                        max_v = max( [max(max(max( v ))), 1e-15] );
                        v_error = max(max(max(abs(v-v_unit)/max_speed)));
                        max_w = max( [max(max(max( abs(w) ))), 1e-15] );
                        w_error = max( [max(max(max(abs(w-w_unit)/max_w))), 1e-15] );
                        max_zeta = max( [max(max(max( zeta ))), 1e-15] );
                        zeta_error = max( [max(max(max(abs(zeta-zeta_unit)/max_zeta))), 1e-15] );
                        max_rho = max( [max(max(max( rho_prime ))), 1e-15] );
                        rho_error = max( [max(max(max(abs(rho_prime-rho_prime_unit)/max_rho))), 1e-15] );
                        
                        max_error = max([round((log10(u_error)))  round((log10(v_error))) round((log10(w_error))) round((log10(zeta_error))) round((log10(rho_error)))],[],'includenan');
                        
                        totalTests = totalTests + 1;
                        if isnan(max_error) || max_error > -10
                            totalErrors = totalErrors + 1;
                            if thesign > 0
                                fprintf('\nFound at large error at lat=%f with %s at +(k,l,j)=(%d,%d,%d):\n', latitude(iLat),apiName,k_loop,l_loop,j0);
                            else
                                fprintf('\nFound at large error at lat=%f  with %s at -(k,l,j)=(%d,%d,%d):\n', latitude(iLat),apiName,k_loop,l_loop,j0);
                            end
                            fprintf('The model solution for (u,v) matches the analytical solution to 1 part in (10^%d, 10^%d) at time t=%d\n', round((log10(u_error))), round((log10(v_error))),t);
                            fprintf('The model solution for (w,zeta,rho_prime) matches the analytical solution to 1 part in (10^%d, 10^%d, 10^%d) at time t=%d\n', round((log10(w_error))), round((log10(zeta_error))), round((log10(rho_error))),t);
                        end
                    end % API loop
                end % sign loop
            end % j loop
        end % l loop
        fprintf('\n');
    end % k loop
end % latitude loop

fprintf('\n***************************************************\n');
if totalErrors > 0
    fprintf('FAILED %d of %d tests.\n',totalErrors, totalTests);
else
    fprintf('PASSED all %d tests.\n', totalTests);
end
