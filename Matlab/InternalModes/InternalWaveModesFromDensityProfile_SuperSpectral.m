%%InternalWaveModesFromDensityProfile Internal wave modes
%
% Given a density profile, this function returns the internal modes and
% their eigenvalues.
%
% 	[F, G, h, N2] = InternalWaveModesFromDensityProfile( rho, z, k,
% 	latitude, norm) returns the internal modes for u & v (F), the
% 	internal modes for w (G), the eigendepths (h), and the buoyancy
% 	frequency used in the calculation (N2). As input you must provide
% 	the density (rho), the depth (z), the wavenumber of wave to be
% 	considered (k), and the norm to be used.
% 
% 	[F, G, h, N2] = InternalWaveModesFromDensityProfile( rho, z, k,
% 	latitude, norm, upper_boundary ) is the same as above, but also
% 	allows you to specify the upper boundary condition. The options are,
%   'rigid_lid' and 'free_surface'.
%
% 	[F, G, h, N2] = InternalWaveModesFromDensityProfile( rho, z, k,
% 	latitude, norm, upper_boundary, fixed_attr ) is the same as above, but
% 	also allows you to fix omega instead of k. The options are, 'fixed_k'
%   or 'fixed_omega'.
% 
% 	The density, rho, and depth, z, can be given as either a row
% 	vector or column vector, but they must agree with each other.
% 
% 	The wavenumber k is 2*pi/wavelength. Typically k=0 for the
% 	geostrophic modes.
% 
% 	The latitude is given in degrees.
% 
% 	The matrices F and G are of size length(z) x length(z). If z is
% 	given as a column vector (with length(z) rows), then each column
% 	of the matrix represents a mode, and vice-a-versa for a row
% 	vector. The G matrix contains the ordered internal modes that for
% 	w. They are the solution to the generalized eigenvalue problem
% 	G_zz = - (1/gh) N^2 G.
% 
% 	The matrices F and G are normalized by given norm, either 1)
% 	max_u, 2) max_w or 3) total_energy. The max_u norm sets the
% 	maximum value of each u,v mode to 1, while the max_w normal sets
% 	the maximum value of each w mode to 1. The total_energy norm uses
% 	\frac{1}{g}\int (N^2 - f_0^2) G^2 \, dz = 1, generally useful
% 	when specify the energy density.
% 
% 	The eigendepths h have units of meters and there is one for each
% 	mode.
% 
% 	The quantity N^2 is the buoyancy frequency.

function [F, G, h, N2] = InternalWaveModesFromDensityProfile_SuperSpectral( rho, z_in, z_out, k, latitude, norm, upper_boundary, fixed_attr )
	
	if strcmp(rho, '--t')
		internal_mode_test,return
    end
    
    if strcmp(rho, '--t2')
		internal_mode_exponential_test,return
	end
	
    if nargin < 8
		fixed_attr = 'fixed_k';
    end
    
	if nargin < 7
		upper_boundary = 'rigid_lid';
	end
		
    if strcmp(fixed_attr, 'fixed_omega')
        omega=k;
    end
    
	if  (~strcmp(upper_boundary, 'free_surface') && ~strcmp(upper_boundary, 'rigid_lid') )
		disp('Invalid upper boundary condition!') 
		return 
	end
	
    if  (~strcmp(norm, 'max_u') && ~strcmp(norm, 'max_w') && ~strcmp(norm, 'total_energy') && ~strcmp(norm, 'const_G_norm') && ~strcmp(norm, 'const_F_norm'))
		disp('Invalid norm!') 
		return 
    end

    g = 9.81;
	f0 = 2*(7.2921e-5)*sin(latitude*pi/180);
    rho_0 = min(rho);
    
    if (z_in(2) - z_in(1)) > 0
        z_in = flip(z_in);
        rho = flip(rho);
        didFlipInputGrid = 1;
    else
        didFlipInputGrid = 0;
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Compute the buoyancy
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
    L_z = max(z_in) - min(z_in);
    z_in_norm_grid = ChebyshevPolynomialsOnGrid( z_in ); % z_in, normalized to [-1, 1] (-1 bottom, 1 top)
    N_points = 4*2*round(pi/(min(abs(diff(z_in_norm_grid)))) + 1);
    z_cheb_grid = cos(((0:N_points-1)')*pi/(N_points-1)); % z, on a chebyshev grid
    % Spline interpolation is necessary, otherwise diff2 exposes piecewise
    % constant bit
    rho_z_norm = interp1(z_in_norm_grid, rho, z_cheb_grid, 'spline'); % rho, interpolated to that grid
    rho_cheb = fct(rho_z_norm);
    rho_cheb(length(z_in):end) = 0;
    
    Diff1 = ChebyshevDifferentiationMatrix( N_points );
    Diff1 = (2/L_z)*Diff1;
    N2_cheb= -(g/rho_0)*Diff1*rho_cheb;
    N2_z_cheb = Diff1*N2_cheb; % end points don't look right.
    
    tmp = ifct(N2_cheb);
    N2_squared_cheb = fct(tmp.*tmp);
    N2_squared_cheb(N_points/2:end) = 0;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Create the stretched grid
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
    s = -g*rho/rho_0;
    L_s = max(s)-min(s);
    N_points = 2*length(z_in); % More is better, but this should be sufficient.
    xi=(0:N_points-1)';
    s_grid = abs(L_s/2)*(cos(xi*pi/(N_points-1))+1) + min(s); % s Chebyshev grid on the s-coordinate
    z_stretched = interp1(s, z_in, s_grid, 'spline'); % z, evaluated on that s grid
    % as defined here, s increases, as z increases
    
    
    N_polys = N_points;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Project N2 & N2_z onto the stretched grid
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
    z_stretched_norm = (2/L_z)*(z_stretched-min(z_in)) - 1; % monotonically decreasing
    t = acos(z_stretched_norm);
    T = zeros(length(z_stretched_norm),N_polys);
    for iPoly=0:(N_polys-1)
       T(:,iPoly+1) = cos(iPoly*t);
    end
    
    N2 = T*N2_cheb(1:N_polys);
    N2_z = T*N2_z_cheb(1:N_polys);
    N2_squared = T*N2_squared_cheb(1:N_polys);
    
    
    %%%%%%%%%%%%%%%%%%%%%
    %
    % Degenerate cases
    %
    %%%%%%%%%%%%%%%%%%%%%
    
    if strcmp(fixed_attr, 'fixed_omega')
        if omega == f0
           k = 0;
           fixed_attr = 'fixed_k';
        elseif omega < f0 % think about this carefully. Really?
            k = 0;
            f0 = 0;
            fixed_attr = 'fixed_k';
        end
    end
    
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Now solve the eigenvalue problem on a Chebyshev Endpoint Grid
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    [T,T_z,T_zz] = ChebyshevPolynomialsOnGrid( s_grid, N_polys );
    
    N0 = 5.23e-3;
    L = 1300;
    N2_b = N0*N0+(2/L)*(s_grid+9.81);
    N2_z_b = (2/L)*(N0*N0+(2/L)*(s_grid+9.81));
    
    if strcmp(fixed_attr, 'fixed_k')
        % The eigenvalue equation is,
        % G_{zz} - K^2 G = \frac{f_0^2 -N^2}{gh_j}G
        % A = \frac{g}{f_0^2 -N^2} \left( \partial_{zz} - K^2*I \right)
        % B = I
           
        A = diag(N2.*N2)*T_zz + diag(N2_z)*T_z - k*k*T;
        B = diag( (f0*f0 - N2)/g )*T;
        
    elseif strcmp(fixed_attr, 'fixed_omega')
        % \frac{f_0^2 - \omega^2}{N^2 - \omega^2} G = k^2 G
        A = diag( (f0*f0 - omega*omega)./(N2 - omega*omega) ) * T_zz;
        B = T;
    else
        disp('Invalid fixed_attr!') 
		return 
    end
    
    % Lower boundary is rigid, G=0
    A(N_polys,:) = T(N_polys,:);
    B(N_polys,:) = 0;
	
	if strcmp(upper_boundary, 'free_surface')
        if strcmp(fixed_attr, 'fixed_k')
        % G_z = \frac{1}{h_j} G at the surface
            A(1,:) = T_z(1,:);
        elseif strcmp(fixed_attr, 'fixed_omega')
           prefactor = (omega*omega - f0*f0)/g;
           A(1,:) = prefactor*T_z(1,:);
        end
		B(1,:) = T(1,:);
	elseif strcmp(upper_boundary, 'rigid_lid')
		A(1,:) = T(1,:);
		B(1,:) = 0;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Now compute the eigenvalues
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
    tic
	[V,D] = eig( A, B );
	toc
    
	[lambda, permutation] = sort(abs(diag(D)),1,'ascend');
    G_in=V(:,permutation);
    
    % Zero pad
%     G_in(N_polys/2:end,:) = 0;
    
    if strcmp(fixed_attr, 'fixed_k')
        h = (1.0 ./ lambda).';
    elseif strcmp(fixed_attr, 'fixed_omega')
        h = ((omega*omega - f0*f0)./(g*lambda)).';
    end
    
    % Create our 'integration' vector for Chebyshev coefficients
%     n = (0:(N_polys-1))';
%     Int1 = -(1+(-1).^n)./(n.*n-1);
%     Int1(2) = 0;
%     Int1 = L/2*Int1;
    
    % This is WRONG, we can't use z_out directly
    %[T] = ChebyshevPolynomialsOnGrid( z_out, N_polys );
    

    
    % Need to use z_in min and max
%     if (z_out(2)-z_out(1) < 0 ) % monotonically decreasing
%         z_out = flip(z_out,1);
%         didFlipOutputGrid = 1;
%     else
%         didFlipOutputGrid = 0; % mono-increasing
%     end

    % Compute N2 on the output grid
    x_norm = (2/(max(z_in)-min(z_in)))*(z_out-min(z_in)) - 1; % mono-increasing
    t = acos(x_norm);
    T = zeros(length(z_out),N_polys);
    for iPoly=0:(N_polys-1)
       T(:,iPoly+1) = cos(iPoly*t);
    end
    N2 = T*N2_cheb(1:N_polys);


    % where in s-coord do we need to evaluate for output points?
    s_out = interp1(z_in, s, z_out)'; % z, evaluated on that s grid
    
    x_norm = (2/L_s)*(s_out-min(s)) - 1; % mono-increasing
    t = acos(x_norm);
    T = zeros(length(z_out),N_polys);
    for iPoly=0:(N_polys-1)
       T(:,iPoly+1) = cos(iPoly*t);
    end
    
    thesign = 1;
%     if (didFlipOutputGrid == 1)
%        T = flip(T,1);
%        z_out = flip(z_out,1); % must be flipped back, because it's used for the norm below.
%        thesign = -1;
%     end
    
    Diff1 = ChebyshevDifferentiationMatrix( N_polys );
    Diff1 = (2/L_s)*Diff1;
    
    F = zeros(length(z_out),N_polys);
    G = zeros(length(z_out),N_polys);
    
    % This alternative method of normalization uses the exact integral of
    % Chebyshev polynomials, rather than the trapz approximation
% 	for j=1:N_polys        
%         if strcmp(norm, 'max_u')
% 			A = max( abs(ifct(h(j) * Diff1 * G_in(:,j))) );
%         elseif strcmp(norm, 'max_w')
% 			A = max( abs(ifct(G_in(:,j))) );
%         elseif strcmp(norm, 'total_energy')
%             J = (1/g) * (N2 - f0*f0) .* ifct(G_in(:,j)) .^ 2;    
% 			A = sqrt(abs(sum(Int1.*fct(J))));
%         end        
%         G(:,j) = T*G_in(:,j) / A;
%         F(:,j) = h(j) * T* Diff1 * G_in(:,j) / A;
% 	end
[~,maxIndexZ] = max(z_out);
tic
    
    
	for j=1:N_polys
        G(:,j) = T*G_in(:,j);
		F(:,j) = h(j) * N2 .* T* Diff1 * G_in(:,j);
%         if canUseFastTransform == 1
%             F(:,j) = ifct( F(:,j) );
%             G(:,j) = ifct( G(:,j) );
%         else
%             F(:,j) = T*F(:,j);
%             G(:,j) = T*G(:,j);
%         end
		if strcmp(norm, 'max_u')
			A = max( abs(F(:,j)) );
			G(:,j) = G(:,j) / A;
			F(:,j) = F(:,j) / A;
		elseif strcmp(norm, 'max_w')
			A = max( abs(G(:,j)) );
			G(:,j) = G(:,j) / A;
			F(:,j) = F(:,j) / A;
		elseif strcmp(norm, 'total_energy') || strcmp(norm, 'const_G_norm')
			A = abs(trapz( z_out, (1/g) * (N2 - f0*f0) .* G(:,j) .^ 2));
			G(:,j) = G(:,j) / sqrt(A);
			F(:,j) = F(:,j) / sqrt(A);
        elseif strcmp(norm, 'const_F_norm')
            A = abs(trapz( z_out, (1/abs(z_in(end)-z_in(1))) .* F(:,j) .^ 2));
            G(:,j) = G(:,j) / sqrt(A);
			F(:,j) = F(:,j) / sqrt(A);
        end
        
        if F(maxIndexZ,j)< 0
            F(:,j) = -F(:,j);
            G(:,j) = -G(:,j);
        end
	end
toc

end

function[]=internal_mode_test
	L_z=-300;
    n = 32;
    xi=(0:n-1)';
    z = abs(L_z/2)*(cos(xi*pi/(n-1))+1) + L_z;
    
	z=(L_z:4.6875:0)';
    z=linspace(L_z,0,n)';
    
    N2 = 1.69e-4;
	rho_0=1025;
    g = 9.81;
	rho = -(N2*rho_0/g)*z + rho_0;
	latitude = 33;
	k=0.1;
    k=0.0;
	
	%[F, G, h, N2] = InternalWaveModesFromDensityProfile_Spectral( rho, z, z, k, latitude, 'max_u', 'free_surface' );
	[F, G, h, N2] = InternalWaveModesFromDensityProfile_SuperSpectral( rho, z, z, k, latitude, 'total_energy', 'rigid_lid' );

    
	figure
	subplot(1,3,1)
	plot(F(:,1:4),z, 'LineWidth', 2)
	ylabel('depth (meters)');
	xlabel('(u,v)-modes');
	
	b = subplot(1,3,2);
	plot(G(:,1:4),z, 'LineWidth', 2)
    title(b, sprintf('Internal Modes with Free Surface\n h = (%.2g, %.2g, %.2g, %.2g)', h(1) , h(2), h(3), h(4) ));
	xlabel('w-modes');
    
    subplot(1,3,3)
	plot(N2,z, 'LineWidth', 2)
    xlim([0.0 1.1*max(N2)])
	xlabel('buoyancy frequency');
end

function[]=internal_mode_exponential_test
	L_z=-5000;
    n = 64;
    xi=(0:n-1)';
    z = abs(L_z/2)*(cos(xi*pi/(n-1))+1) + L_z;
    
	z=(L_z:4.6875:0)';
    z=linspace(L_z,0,n)';
    
    
    N0 = 5.23e-3;
    b = 1300;
    g = 9.81;
    rho0=1025;
    rho = rho0*N0*N0*b/(2*g)*(1-exp(2*z/1300))+rho0;
	latitude = 33;
	k=0.1;
    k=0.0;
	
	%[F, G, h, N2] = InternalWaveModesFromDensityProfile_Spectral( rho, z, z, k, latitude, 'max_u', 'free_surface' );
	[F, G, h, N2] = InternalWaveModesFromDensityProfile_SuperSpectral( rho, z, z, k, latitude, 'max_u', 'rigid_lid' );

    
	figure
	subplot(1,3,1)
	plot(F(:,1:4),z, 'LineWidth', 2)
	ylabel('depth (meters)');
	xlabel('(u,v)-modes');
	
	b = subplot(1,3,2);
	plot(G(:,1:4),z, 'LineWidth', 2)
    title(b, sprintf('Internal Modes with Free Surface\n h = (%.2g, %.2g, %.2g, %.2g)', h(1) , h(2), h(3), h(4) ));
	xlabel('w-modes');
    
    subplot(1,3,3)
	plot(N2,z, 'LineWidth', 2)
    xlim([0.0 1.1*max(N2)])
	xlabel('buoyancy frequency');
end