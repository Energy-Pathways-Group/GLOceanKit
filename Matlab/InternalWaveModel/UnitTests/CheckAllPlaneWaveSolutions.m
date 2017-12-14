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
Nz = 2;

latitude = 31;
N0 = 5.2e-3/2; % Choose your stratification 7.6001e-04

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if 1 == 1
    wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);
    rho0 = wavemodel.rho0;
    g = 9.81;
else
    rho0 = 1025; g = 9.81;
    rho = @(z) -(N0*N0*rho0/g)*z + rho0;
    z = (Lz/Nz)*(0:Nz-1)' - Lz;
    wavemodel = InternalWaveModelArbitraryStratification([Lx, Ly, Lz], [Nx, Ny, Nz], rho, z, Nz, latitude, 'method','wkbSpectral','nEVP',128);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create a single plane-wave with the model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

U = 0.01; % m/s
phi = 0*0.232; % random, just to see if its doing the right thing

totalErrors = 0;
totalTests = 0;
for k_loop=(-Nx/2 + 1):1:(Nx/2-1)
    for l_loop=(-Ny/2 + 1):1:(Ny/2-1)
        fprintf('(k0,l0)=(%d,%d) ',k_loop,l_loop);
        for j0=1:Nz/2
            for sign=-1:2:1
                
                f0 = wavemodel.f0;
                k = wavemodel.k;
                l = wavemodel.l;
                h = wavemodel.h;
                x = wavemodel.x;
                y = wavemodel.y;
                X = wavemodel.X;
                Y = wavemodel.Y;
                Z = wavemodel.Z;
                
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
                
                if 1 == 1
                    period = wavemodel.InitializeWithPlaneWave(k_loop,l_loop,j0,U,sign);
                    omega = sign*2*pi/period;
                    kk = k(k0+1);
                    ll = l(l0+1);
                else
                    [omega,kk,ll] = wavemodel.SetGriddedWavesWithWavemodes(k_loop,l_loop,j0,phi,U,sign);
                    period = 2*pi/abs(omega);
                end
                
                alpha=atan2(ll,kk);
                K = sqrt( kk^2 + ll^2);
                
                t = 4*86400;
                [u,v,w,zeta] = wavemodel.VariableFieldsAtTime(t,'u','v','w','zeta');
                rho = wavemodel.DensityFieldAtTime(t);
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %
                % Create a single plane-wave with the known analytical solution
                %
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%              

                if latitude == 0
                    f0OverOmega = 0;
                    if K == 0
                        kOverOmega = 0;
                    else
                        kOverOmega = K/omega;
                    end
                else
                    f0OverOmega = f0/omega;
                    kOverOmega = K/omega;
                end
                theta = kk*X + ll*Y + omega*t + phi;
                u_unit = U*(cos(alpha)*cos( theta ) + f0OverOmega*sin(alpha)*sin( theta )).*cos(m*Z);
                v_unit = U*(sin(alpha)*cos( theta ) - f0OverOmega*cos(alpha)*sin( theta )).*cos(m*Z);
                w_unit = (U*K/m) * sin( theta ) .* sin(m*Z);
                zeta_unit = -(U/m) * kOverOmega * cos( theta ) .* sin(m*Z);
                rho_unit = rho0 -(N0*N0*rho0/g)*reshape(wavemodel.z,1,1,[]) -(rho0/g)*N0*N0*(U*K/m/omega) * cos( theta ) .* sin(m*Z);
                
                if any(any(any(isnan(zeta_unit))))
                   error('nan') 
                end
                
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
                max_rho = max( [max(max(max( rho ))), 1e-15] );
                rho_error = max( [max(max(max(abs(rho-rho_unit)/max_rho))), 1e-15] );
                
                max_error = max([round((log10(u_error)))  round((log10(v_error))) round((log10(w_error))) round((log10(zeta_error))) round((log10(rho_error)))],[],'includenan');
                
                totalTests = totalTests + 1;
                if isnan(max_error) || max_error > -3
                    totalErrors = totalErrors + 1;
                    if sign > 0
                        fprintf('\nFound at large error at +(k,l,j)=(%d,%d,%d):\n',k_loop,l_loop,j0);
                    else
                        fprintf('\nFound at large error at -(k,l,j)=(%d,%d,%d):\n',k_loop,l_loop,j0);
                    end
                    fprintf('The model solution for (u,v) matches the analytical solution to 1 part in (10^%d, 10^%d) at time t=%d\n', round((log10(u_error))), round((log10(v_error))),t);
                    fprintf('The model solution for (w,zeta,rho) matches the analytical solution to 1 part in (10^%d, 10^%d, 10^%d) at time t=%d\n', round((log10(w_error))), round((log10(zeta_error))), round((log10(rho_error))),t);
                end
            end
        end
    end
    fprintf('\n');
end

fprintf('\n***************************************************\n');
if totalErrors > 0
    fprintf('FAILED %d of %d tests.\n',totalErrors, totalTests);
else
    fprintf('PASSED all %d tests.\n', totalTests);
end
