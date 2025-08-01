% This test suite tests the following methods:
% -diffX
% -diffY
% -diffZF
% -diffZG

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   Fourier derivatives
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
fprintf('Testing Fourier transform:\n')

Nx = 16;
Ny = 8;
Nz = 16+1;

Lx = 1.0;
Ly = 10.0;
Lz = 4.0;

latitude = 25;
N0 = 5.2e-3; % Choose your stratification 7.6001e-04
wvt = WVTransformConstantStratification([Lx, Ly, Lz], [Nx, Ny, Nz], N0, latitude=latitude);

nu = 2*pi/Lx;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nth-derivative, cosine, 1st-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y,Z] = wvt.xyzGrid;
f = @(omega) cos(omega*X) .* cos(nu*Y);

% range of frequencies to test
iOmega = 2*pi*(1:floor(Nx/2))'/Lx;

for iDerivative=1:4
    if iDerivative==1
        Df_analytical = @(omega) -omega*sin(omega*X).*cos(nu*Y);
    elseif iDerivative==2
        Df_analytical = @(omega) -(omega^2)*cos(omega*X).*cos(nu*Y);
    elseif iDerivative==3
        Df_analytical = @(omega) (omega^3)*sin(omega*X).*cos(nu*Y);
    elseif iDerivative==4
        Df_analytical = @(omega) (omega^4)*cos(omega*X).*cos(nu*Y);
    end

    Df_numerical = @(u) wvt.diffX(u,iDerivative);
    testname = sprintf('Fourier differentiation of cosine (numDerivs=%d), 1st dimension',iDerivative);
    ReportErrors(f,Df_analytical,Df_numerical,testname,iOmega);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nth-derivative, cosine, 2nd-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y,Z] = wvt.xyzGrid;
f = @(omega) cos(nu*X) .* cos(omega*Y);

% range of frequencies to test
iOmega = 2*pi*(1:floor(Ny/2))'/Ly;

for iDerivative=1:4
    if iDerivative==1
        Df_analytical = @(omega) -omega*sin(omega*Y).*cos(nu*X);
    elseif iDerivative==2
        Df_analytical = @(omega) -(omega^2)*cos(omega*Y).*cos(nu*X);
    elseif iDerivative==3
        Df_analytical = @(omega) (omega^3)*sin(omega*Y).*cos(nu*X);
    elseif iDerivative==4
        Df_analytical = @(omega) (omega^4)*cos(omega*Y).*cos(nu*X);
    end

    Df_numerical = @(u) wvt.diffY(u,iDerivative);
    testname = sprintf('Fourier differentiation of cosine (numDerivs=%d), 2nd dimension',iDerivative);
    ReportErrors(f,Df_analytical,Df_numerical,testname,iOmega);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nth cosine derivative, 3rd-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y,Z] = wvt.xyzGrid;
f = @(omega) cos(nu*X).*cos(omega*Z);

% range of frequencies to test
df = 1/((Nz-1)*(wvt.z(2)-wvt.z(1)));
iOmega = pi*df*(0:(Nz-1))';

for iDerivative=1:4
    if iDerivative==1
        Df_analytical = @(omega) -omega*sin(omega*Z).*cos(nu*X);
    elseif iDerivative==2
        Df_analytical = @(omega) -(omega^2)*cos(omega*Z).*cos(nu*X);
    elseif iDerivative==3
        Df_analytical = @(omega) (omega^3)*sin(omega*Z).*cos(nu*X);
    elseif iDerivative==4
        Df_analytical = @(omega) (omega^4)*cos(omega*Z).*cos(nu*X);
    end

    Df_numerical = @(u) wvt.diffZF(u,iDerivative);
    testname = sprintf('F (cosine) differentiation of cosine (numDerivs=%d), 3rd dimension',iDerivative);
    ReportErrors(f,Df_analytical,Df_numerical,testname,iOmega);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% nth sine derivative, 3rd-dimension
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[X,Y,Z] = wvt.xyzGrid;
f = @(omega) cos(nu*X).*sin(omega*Z);

% range of frequencies to test
df = 1/((Nz-1)*(wvt.z(2)-wvt.z(1)));
iOmega = pi*df*(0:(Nz-1))';

for iDerivative=1:4
    if iDerivative==1
        Df_analytical = @(omega) omega*cos(omega*Z).*cos(nu*X);
    elseif iDerivative==2
        Df_analytical = @(omega) -(omega^2)*sin(omega*Z).*cos(nu*X);
    elseif iDerivative==3
        Df_analytical = @(omega) -(omega^3)*cos(omega*Z).*cos(nu*X);
    elseif iDerivative==4
        Df_analytical = @(omega) (omega^4)*sin(omega*Z).*cos(nu*X);
    end

    Df_numerical = @(u) wvt.diffZG(u,iDerivative);
    testname = sprintf('G (sin) differentiation of cosine (numDerivs=%d), 3rd dimension',iDerivative);
    ReportErrors(f,Df_analytical,Df_numerical,testname,iOmega);
end


