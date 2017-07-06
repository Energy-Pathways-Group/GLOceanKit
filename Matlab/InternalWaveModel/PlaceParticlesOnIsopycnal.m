function zIsopycnal = PlaceFloatsOnIsopycnal(x,y,z,X,Y,Z,Rho, drhobardz, interpolationMethod,tolerance)
% (x,y,z) are the positions of the floats. Any floats with the same value
% of z will be moved to the same isopycnal.
%
% (X,Y,Z) are the domain dimensions (in ndgrid format) for Rho. By
% assumption the domain is periodic in (x,y).
%
% interpolation method should be either 'spline' or 'linear', try 'spline'.
%
% tolerance is in meters, 1e-8 should be fine.

% There must be a better way! What a horrible hack.
if (X(2,1,1) - X(1,1,1) > 0)
    xdim = 1;
    xstride = 1;
elseif (X(1,2,1) - X(1,1,1) > 0)
    xdim = 2;
    xstride = size(X,1);
elseif (X(1,1,2) - X(1,1,1) > 0)
    xdim = 3;
    xstride = size(X,1)*size(X,2);
end

if (Y(2,1,1) - Y(1,1,1) > 0)
    ydim = 1;
    ystride = 1;
elseif (Y(1,2,1) - Y(1,1,1) > 0)
    ydim = 2;
    ystride = size(Y,1);
elseif (Y(1,1,2) - Y(1,1,1) > 0)
    ydim = 3;
    ystride = size(Y,1)*size(Y,2);
end

if (Z(2,1,1) - Z(1,1,1) > 0)
    zdim = 1;
    zstride = 1;
elseif (Z(1,2,1) - Z(1,1,1) > 0)
    zdim = 2;
    zstride = size(Y,1);
elseif (Z(1,1,2) - Z(1,1,1) > 0)
    zdim = 3;
    zstride = size(Y,1)*size(Y,2);
end

Nx = size(X,xdim);
Ny = size(Y,ydim);
Nz = size(Z,zdim);
dx = X(2*xstride) - X(1*xstride);
dy = Y(2*ystride) - Y(1*ystride);
Lx = dx*(Nx+1); % The "plus one" is because it's periodic.
Ly = dy*(Ny+1);

zCoordinate = Z((1:Nz)*zstride);
drhobardzFunction = @(z) interp1(zCoordinate,drhobardz,z,'spline');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Wrap the particle positions, as necessary
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% (x,y) are periodic for the gridded solution
x_tilde = mod(x,Lx);
y_tilde = mod(y,Ly);

x_index = floor(x_tilde/dx);
y_index = floor(y_tilde/dy);

% Identify the particles along the interpolation boundary
if strcmp(interpolationMethod,'spline')
    S = 3+1; % cubic spline, plus buffer
elseif strcmp(interpolationMethod,'linear')
    S = 1+1;
end
bp = x_index < S-1 | x_index > Nx-S | y_index < S-1 | y_index > Ny-S;

% then do the same for particles that fall along the boundary.
x_tildeS = mod(x(bp)+S*dx,Lx);
y_tildeS = mod(y(bp)+S*dy,Ly);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Create the interpolants
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RhoInterp = griddedInterpolant(X,Y,Z, Rho,interpolationMethod);
if any(bp)
    shiftVector = zeros(1,3);
    shiftVector(xdim) = S;
    shiftVector(ydim) = S;
    RhoInterpS = griddedInterpolant(X,Y,Z, circshift(Rho,shiftVector),interpolationMethod);
end

% Now let's place the floats along an isopycnal.
zIsopycnal = z;

zUnique = unique(z);
rho = zeros(size(zIsopycnal));
for zLevel = 1:length(zUnique)
    iterations = 0;
    zLevelIndices = (z==zUnique(zLevel));
    
    nonboundaryIndices = zLevelIndices & ~bp;
    rho(nonboundaryIndices) = RhoInterp(x_tilde(nonboundaryIndices),y_tilde(nonboundaryIndices),zIsopycnal(nonboundaryIndices));
    if any(bp)
        boundaryIndices = zLevelIndices & bp;
        rho(boundaryIndices) = RhoInterpS(x_tildeS(boundaryIndices),y_tildeS(boundaryIndices),zIsopycnal(boundaryIndices));
    end
        
    dRho = rho(zLevelIndices) - mean(rho(zLevelIndices));
    dz = -dRho./drhobardzFunction(zIsopycnal(zLevelIndices));
    zIsopycnal(zLevelIndices) = zIsopycnal(zLevelIndices)+dz;
    
    while( max(abs(dz)) > tolerance && iterations < 20 )
        rho(nonboundaryIndices) = RhoInterp(x_tilde(nonboundaryIndices),y_tilde(nonboundaryIndices),zIsopycnal(nonboundaryIndices));
        if any(bp)
            rho(boundaryIndices) = RhoInterpS(x_tildeS(boundaryIndices),y_tildeS(boundaryIndices),zIsopycnal(boundaryIndices));
        end
        
        dRho = rho(zLevelIndices) - mean(rho(zLevelIndices));
        dz = -dRho./drhobardzFunction(zIsopycnal(zLevelIndices));
        zIsopycnal(zLevelIndices) = zIsopycnal(zLevelIndices)+dz;
        
        iterations = iterations + 1;
    end
    
    if (iterations == 20)
        fprintf('Exceed maximum number of iterations\n');
    end
    
    fprintf('All floats are within %.2g meters of the isopycnal at z=%.1f meters\n',max(abs(dz)),mean(z(zLevelIndices)) )
end
end