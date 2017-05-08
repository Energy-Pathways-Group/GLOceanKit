Lx = 15e3;
Ly = 15e3;
Lz = 5000;

Nx = 8;
Ny = 8;
Nz = 5; % Must include end point to advect at the surface, so use 2^N + 1

latitude = 31;
N0 = 5.2e-3; % Choose your stratification

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Initialize the wave model
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


wavemodel = InternalWaveModelConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], latitude, N0);

% Our measure of error.
error = @(u,u_unit) max( [max(max(max(abs(u-u_unit)/(max(max(max( u_unit ))) - min(min(min( u_unit )))) ))), 1e-15]);

% We test the interpolation method by checking to see that the velocity
% field on the grid points, lines up with the FFT solution at the same
% point. Although it doesn't prove that the spectral interpolation is
% working off-grid, it is likely.
x_float = reshape(wavemodel.X,[],1);
y_float = reshape(wavemodel.Y,[],1);
z_float = reshape(wavemodel.Z,[],1);
p0 = [x_float y_float z_float];

for k_loop=(-Nx/2 + 1):1:(Nx/2-1)
    for l_loop=(-Ny/2 + 1):1:(Ny/2-1)
        fprintf('(k0,l0)=(%d,%d) ',k_loop,l_loop);
        for j0=1:Nz/2
            for sign=-1:2:1
                U = 0.01;
                period = wavemodel.InitializeWithPlaneWave(k_loop,l_loop,j0,U,sign);

                [U,V,W] = wavemodel.VelocityFieldAtTime(0);
                [~,ZETA] = wavemodel.VerticalFieldsAtTime(0);
                u = reshape(U,[],1);
                v = reshape(V,[],1);
                w = reshape(W,[],1);
                zeta = reshape(ZETA,[],1);
                
                [u_spec, v_spec, w_spec] = wavemodel.VelocityAtTimePosition(0, p0(:,1), p0(:,2), p0(:,3), 'exact');
                [zeta_spec] = wavemodel.ZetaAtTimePosition(0, p0(:,1), p0(:,2), p0(:,3));
                
                
                
                u_error = error(u_spec,u);
                v_error = error(v_spec,v);
                w_error = error(w_spec,w);
                zeta_error = error(zeta_spec,zeta);
                
                max_error = max([round((log10(u_error)))  round((log10(v_error))) round((log10(w_error))) round((log10(zeta_error)))]);
                
                if max_error > -3
                    if sign > 0
                        fprintf('\nFound at large error at +(k,l,j)=(%d,%d,%d):\n',k_loop,l_loop,j0);
                    else
                        fprintf('\nFound at large error at -(k,l,j)=(%d,%d,%d):\n',k_loop,l_loop,j0);
                    end
                    fprintf('The spectral interpolation for (u,v,w,zeta) matches FFT solution to 1 part in (10^%d, 10^%d, 10^%d, 10^%d)\n', round((log10(u_error))), round((log10(v_error))), round((log10(w_error))), round((log10(zeta_error))) );
                end
            end
        end
    end
end
fprintf('\n');

wavemodel.FillOutWaveSpectrum();
wavemodel.InitializeWithGMSpectrum(1.0);

[U,V,W] = wavemodel.VelocityFieldAtTime(0);
[~,ZETA] = wavemodel.VerticalFieldsAtTime(0);
u = reshape(U,[],1);
v = reshape(V,[],1);
w = reshape(W,[],1);
zeta = reshape(ZETA,[],1);

[u_spec, v_spec, w_spec] = wavemodel.VelocityAtTimePosition(0, p0(:,1), p0(:,2), p0(:,3), 'exact');
[zeta_spec] = wavemodel.ZetaAtTimePosition(0, p0(:,1), p0(:,2), p0(:,3));
u_error = error(u_spec,u);
v_error = error(v_spec,v);
w_error = error(w_spec,w);
zeta_error = error(zeta_spec,zeta);

max_error = max([round((log10(u_error)))  round((log10(v_error))) round((log10(w_error))) round((log10(zeta_error)))]);
fprintf('The spectral interpolation for (u,v,w,zeta) matches FFT solution to 1 part in (10^%d, 10^%d, 10^%d, 10^%d)\n', round((log10(u_error))), round((log10(v_error))), round((log10(w_error))), round((log10(zeta_error))) );
