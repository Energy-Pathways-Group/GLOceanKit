function [F, G, h, N2] = InternalWaveModesFromDensityProfile_FiniteDifference( rho, z, k, latitude, norm, upper_boundary, orderOfAccuracy, fixed_attr )
%% InternalWaveModesFromDensityProfile_FiniteDifference
%
% Given a density profile, this function returns the internal modes and
% their eigenvalues using finite differencing.
%
% 	[F, G, h, N2] = InternalWaveModesFromDensityProfile( rho, z, k,
% 	latitude, norm) returns the internal modes for u & v (F), the internal
% 	modes for w (G), the eigendepths (h), and the buoyancy frequency used
% 	in the calculation (N2). As input you must provide the density (rho),
% 	the depth (z), the wavenumber of wave to be considered (k), and the
% 	norm to be used.
% 
% 	[F, G, h, N2] = InternalWaveModesFromDensityProfile( rho, z, k,
% 	latitude, norm, upper_boundary ) is the same as above, but also allows
% 	you to specify the upper boundary condition. The options are,
%   'rigid_lid' and 'free_surface'.
%
% 	[F, G, h, N2] = InternalWaveModesFromDensityProfile( rho, z, k,
% 	latitude, norm, upper_boundary, orderOfAccuracy ) is the same as above,
%   but also allows you to specify the order of the finite difference
%   matrix. Minimum of 2, default is currently 4.
%
% 	[F, G, h, N2] = InternalWaveModesFromDensityProfile( rho, z, k,
% 	latitude, norm, upper_boundary, orderOfAccuracy, fixed_attr ) is the
% 	same as above, but also allows you to fix omega instead of k. The
% 	options are, 'fixed_k' or 'fixed_omega'.
% 
% 	The density, rho, and depth, z, can be given as either a row vector or
% 	column vector, but they must agree with each other.
% 
% 	The wavenumber k is 2*pi/wavelength. Typically k=0 for the geostrophic
% 	modes.
% 
% 	The latitude is given in degrees.
% 
% 	The matrices F and G are of size length(z) x length(z). If z is given
% 	as a column vector (with length(z) rows), then each column of the
% 	matrix represents a mode, and vice-a-versa for a row vector. The G
% 	matrix contains the ordered internal modes that for w. They are the
% 	solution to the generalized eigenvalue problem G_zz = - (1/gh) N^2 G.
% 
% 	The matrices F and G are normalized by given norm, either 1) max_u, 2)
% 	max_w or 3) total_energy. The max_u norm sets the maximum value of each
% 	u,v mode to 1, while the max_w normal sets the maximum value of each w
% 	mode to 1. The total_energy norm uses \frac{1}{g}\int (N^2 - f_0^2) G^2
% 	\, dz = 1, generally useful when specify the energy density.
% 
% 	The eigendepths h have units of meters and there is one for each mode.
% 
% 	The quantity N^2 is the buoyancy frequency.
	
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
		orderOfAccuracy = 4;
    end
    
	if nargin < 6
		upper_boundary = 'rigid_lid';
	end
	
	if abs(std(diff(z))/mean(diff(z))) < 1e-6
        %disp('This is an evenly spaced grid.') 
	else
        z_norm = z - min(z);
        L = z_norm(1);
        z_norm = z_norm*(2/L)-1;
        n = length(z);
        xi=(0:n-1)';
        z_compare = cos(xi*pi/(n-1));
        z_diff = z_norm-z_compare;
        if max(abs(z_diff)) < 1e-6 
            disp('This is a Chebyshev endpoint grid, and should probably be solved using InternalWaveModesFromDensityProfile_Spectral. Proceeding anyway.');
        end
	end
	
    if strcmp(fixed_attr, 'fixed_omega')
        omega=k;
    end
    
	if  (~strcmp(upper_boundary, 'free_surface') && ~strcmp(upper_boundary, 'rigid_lid') )
		disp('Invalid upper boundary condition!') 
		return 
	end
	
	if  (~strcmp(norm, 'max_u') && ~strcmp(norm, 'max_w') && ~strcmp(norm, 'total_energy'))
		disp('Invalid norm!') 
		return 
	end

	didFlip = 0;
	if (z(2) - z(1)) < 0
		if (isrow(z))
			z=fliplr(z);
		else
			z=flipud(z);
		end
		if (isrow(rho))
			rho=fliplr(rho);
		else
			rho=flipud(rho);
		end
		didFlip=1;
	end

	N = length(z);
	g = 9.81;
	f0 = 2*(7.2921e-5)*sin(latitude*pi/180);
    
	if orderOfAccuracy < 2
        orderOfAccuracy = 2;
	elseif orderOfAccuracy > N
        orderOfAccuracy = N;
	end
        
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
	% First order differentiation matrix.
	% Neumann boundary conditions.
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	% The operator D is the first order, central difference operator
	Diff1 = FiniteDifferenceMatrix(1, z, 1, 1, orderOfAccuracy);

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Compute the buoyancy
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	%rho_0=trapz( z, rho)/(z(end)-z(1));
    rho_0=min(rho);

	% Convert density to buoyancy. 
	Bouyancy=-g*rho/rho_0; 

	% N^2 is computed assuming that d^2rho/dz^2=0 at the end points.
	N2=Diff1*Bouyancy;

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Second order differentiation matrix with unspecified boundary conditions.
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%

	if strcmp(upper_boundary, 'free_surface')
        rightBCDerivs = 1;
	else
        rightBCDerivs = 0;
	end
	Diff2 = FiniteDifferenceMatrix(2, z, 0, rightBCDerivs, orderOfAccuracy);
	
    if strcmp(fixed_attr, 'fixed_k')
        % The eigenvalue equation is,
        % G_{zz} - K^2 G = \frac{f_0^2 -N^2}{gh_j}G
        % A = \frac{g}{f_0^2 -N^2} \left( \partial_{zz} - K^2*I \right)
        % B = I
        A = diag(-g./(N2 - f0*f0)) * (Diff2 - k*k*eye(N));
        B = eye(N);
    elseif strcmp(fixed_attr, 'fixed_omega')
        % \frac{f_0^2 - \omega^2}{N^2 - \omega^2} G = k^2 G
        A = diag( (f0*f0 - omega*omega)./(N2 - omega*omega) ) * Diff2;
        B = eye(N);
    else
        disp('Invalid fixed_attr!') 
		return 
    end
    
    % G=0 at the bottom
    A(1,:) = Diff2(1,:);
    B(1,:) = 0;
    
    % G=0 or G_z = \frac{1}{h_j} G at the surface, depending on the BC
    A(end,:) = Diff2(end,:);
	if strcmp(upper_boundary, 'free_surface')
        if strcmp(fixed_attr, 'fixed_k')
            % G_z = \frac{1}{h_j} G at the surface
             A(end,:) = Diff2(end,:);
        elseif strcmp(fixed_attr, 'fixed_omega')
            prefactor = (omega*omega - f0*f0)/g;
            A(end,:) = prefactor*Diff2(end,:);
        end
		B(end,end)=1;
	elseif strcmp(upper_boundary, 'rigid_lid')
        % G=0 at the surface (note we chose this BC when creating Diff2)
        A(end,:) = Diff2(end,:);
		B(end,end)=0;
	end

	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	%
	% Now compute the eigenvalues
	%
	%%%%%%%%%%%%%%%%%%%%%%%%%%%
	
	[V,D] = eig( A, B );
	
	[lambda, permutation] = sort(real(diag(D)),1,'ascend');
    G=V(:,permutation);
    
    if strcmp(fixed_attr, 'fixed_k')
        h = (1.0 ./ lambda).';
    elseif strcmp(fixed_attr, 'fixed_omega')
        h = ((omega*omega - f0*f0)./(g*lambda)).';
    end
	
    [~,maxIndexZ] = max(z);
    F = zeros(size(G));
	for j=1:length(G(1,:))
		F(:,j) = h(j) * Diff1 * G(:,j);
		if strcmp(norm, 'max_u')
			A = max( abs(F(:,j)) );
			G(:,j) = G(:,j) / A;
			F(:,j) = F(:,j) / A;
		elseif strcmp(norm, 'max_w')
			A = max( abs(G(:,j)) );
			G(:,j) = G(:,j) / A;
			F(:,j) = F(:,j) / A;
		elseif strcmp(norm, 'total_energy')
			A = trapz( z, (1/g) * (N2 - f0*f0) .* G(:,j) .^ 2);
			G(:,j) = G(:,j) / sqrt(A);
			F(:,j) = F(:,j) / sqrt(A);
        end
        
        if F(maxIndexZ,j)< 0
            F(:,j) = -F(:,j);
            G(:,j) = -G(:,j);
        end
	end


	if (didFlip == 1)
		if (isrow(N2))
			N2=fliplr(N2);
		else
			N2=flipud(N2);
		end
	
		G = flip(G,1);
		F = flip(F,1);
	end
end

function[]=internal_mode_test
	L_z=-300;
	z=(L_z:0.5:0)';
	rho_0=1024.63;
	rho_L=1026.45;
	rho = ((rho_0 - rho_L)/(0-L_z)).*z + rho_0;
	latitude = 31;
	k=0.1;
	
	[F, G, h, N2] = InternalWaveModesFromDensityProfile_FiniteDifference( rho, z, k, latitude, 'max_u', 'free_surface', 8 );

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
%     k=0.0;
	
	[F, G, h, N2] = InternalWaveModesFromDensityProfile_FiniteDifference( rho, z, k, latitude, 'max_u', 'rigid_lid' );

    
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