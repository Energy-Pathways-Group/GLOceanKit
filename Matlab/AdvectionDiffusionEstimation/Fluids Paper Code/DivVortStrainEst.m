function [div, vort, nstrain, sstrain, u_bg, v_bg, E, F] = DivVortStrainEst( x, y, u, v, dt, isDeltaZero, isZetaZero, isSigma_nZero, isSigma_sZero)

% x and y are position data, they should be matrices of size txd, where t
% is the number of time points and d is the number of drifters.

% delta, zeta, signa_n and sigma_s should be entered as either 1 or 0. This
% is a check - is the variable known to be zero? 1 means the variable is
% known to be zero, 0 means the variable is not known to be zero and is to
% be estimated.

%
%
%
%
%First check if delta, zeta, signa_n and sigma_s are all either zero or 1.


if (isDeltaZero ~= 0 && isDeltaZero ~= 1)
    
    msg = 'Error occurred. Please enter either 1 or 0 for delta, where 1 means that the divergence is known to be zero, and 0 means it is to be estimated.';
    error(msg)

elseif (isZetaZero ~= 0 && isZetaZero ~= 1)
        
    msg = 'Error occurred. Please enter either 1 or 0 for zeta, where 1 means that the vorticity is known to be zero, and 0 means it is to be estimated.';
    error(msg)
    
elseif (isSigma_nZero ~= 0 && isSigma_nZero ~= 1)
        
    msg = 'Error occurred. Please enter either 1 or 0 for sigma_n, where 1 means that the normal strain is known to be zero, and 0 means it is to be estimated.';
    error(msg)
    
elseif (isSigma_sZero ~= 0 && isSigma_sZero ~= 1)
        
    msg = 'Error occurred. Please enter either 1 or 0 for sigma_s, where 1 means that the sheer strain is known to be zero, and 0 means it is to be estimated.';
    error(msg)
    
% elseif (isDeltaZero == 1 && isZetaZero ==1 && isSigma_nZero == 1 && isSigma_sZero == 1) % if all variables are zero then the solution is known and estimation isn't required.
% 
% msg = 'Error occurred. All variables are known to be zero. No estimation is required.';
% error(msg)
%     
else
    %
    %
    %
    %
    % Now we can start defining variables to use throughout the function

    %Initialise the variables
    onesvec = ones([length(u(1,:))*length(u(:,1)),1]);
    R= zeros(length(u(1,:))*length(u(:,1)),3);

    U=zeros(length(u(1,:))*length(u(:,1)),1);
    V=zeros(length(u(1,:))*length(u(:,1)),1);


    %Fill the matrices/vectors with data

    R(:,:) = [onesvec, x(:), y(:)];

    U=u(:);
    V=v(:);

    %We now change the above matrix and vectors to allow for setting
    %variables to be zero. We will be estimating A and B simultaneously, so
    %put everything into a big matrix/vector

    R = [R, zeros(size(R)) ; zeros(size(R)), R];

    UV = [U ; V];

    %Edit R for zero variables - The R matrix must be ammended to allow
    %some variables to be known to be zero. This changes for each
    %combination of known zeros, so the possible options are all defined
    %below:
    
    n=size(R,1);
    
    %no variables are known to be zero
    if (isDeltaZero == 0 && isZetaZero == 0 && isSigma_nZero == 0 && isSigma_sZero == 0)
        
        AB = ((R'*R)\R')*UV;
        ux=AB(2);
        uy=AB(3);
        vx=AB(5);
        vy=AB(6);
        u_bg = AB(1); 
    v_bg = AB(4);
    
    %%%%%%%% one variable is zero %%%%%%%% %
    %divergence free only    
    elseif (isDeltaZero == 1 && isZetaZero == 0 && isSigma_nZero == 0 && isSigma_sZero == 0)
        
        R(n/2+1:n,6) = 0;
        R(n/2+1:n,2) = -y(:);
        
        R = R(:,1:5); % remove empty columns
        
        %estimate variables
        AB = ((R'*R)\R')*UV;
        ux=AB(2);
        uy=AB(3);
        vx=AB(5);
        vy = -ux;
        u_bg = AB(1); 
    v_bg = AB(4);
        
    %vorticity free only
    elseif (isDeltaZero == 0 && isZetaZero == 1 && isSigma_nZero == 0 && isSigma_sZero == 0)
        
        R(n/2+1:n,5) = 0;
        R(n/2+1:n,3) = x(:);
        
        R = [R(:,1:4) R(:,6)]; % remove empty columns
        
        %estimate variables
        AB = ((R'*R)\R')*UV;
        ux=AB(2);
        uy=AB(3);
        vx=uy;
        vy=AB(5); % no entry for vx, so vy becomes the 5th entry
        u_bg = AB(1); 
    v_bg = AB(4);
        
    %normal strain = 0 only
    elseif (isDeltaZero == 0 && isZetaZero == 0 && isSigma_nZero == 1 && isSigma_sZero == 0)        
        
        R(n/2+1:n,6) = 0;
        R(n/2+1:n,2) = y(:);
        
        R = R(:,1:5); % remove empty columns
        
        %estimate variables
        AB = ((R'*R)\R')*UV;  
        ux=AB(2);
        uy=AB(3);
        vx=AB(5);
        vy = ux;
        u_bg = AB(1); 
    v_bg = AB(4);
        
    %shear strain = 0 only
    elseif (isDeltaZero == 0 && isZetaZero == 0 && isSigma_nZero == 0 && isSigma_sZero == 1)        
        
        R(n/2+1:n,5) = 0;
        R(n/2+1:n,3) = -x(:);
        
        R = [R(:,1:4) R(:,6)]; % remove empty columns        
        
        %estimate variables
        AB = ((R'*R)\R')*UV;
        ux=AB(2);
        uy=AB(3);
        vx=-uy;
        vy=AB(5); % no entry for vx, so vy becomes the 5th entry
        u_bg = AB(1); 
    v_bg = AB(4);
    
    %%%%%%%% two variables are zero %%%%%%%% %
    %divergence and vorticity = 0 
    elseif (isDeltaZero == 1 && isZetaZero == 1 && isSigma_nZero == 0 && isSigma_sZero == 0)  
        
        R(n/2+1:n,6) = 0;
        R(n/2+1:n,2) = -y(:); 
        R(n/2+1:n,5) = 0;
        R(n/2+1:n,3) = x(:);
        
        R = R(:,1:4); % remove empty columns
        
        %estimate variables
        AB = ((R'*R)\R')*UV;
        ux=AB(2);
        uy=AB(3);
        vx=uy;
        vy=-ux;
        u_bg = AB(1); 
    v_bg = AB(4);
        
    %normal and shear strain = 0 
    elseif (isDeltaZero == 0 && isZetaZero == 0 && isSigma_nZero == 1 && isSigma_sZero == 1)  
        
        R(n/2+1:n,6) = 0;
        R(n/2+1:n,2) = y(:); 
        R(n/2+1:n,5) = 0;
        R(n/2+1:n,3) = -x(:);
        
        R = R(:,1:4); % remove empty columns
        
        %estimate variables
        AB = ((R'*R)\R')*UV;
        ux=AB(2);
        uy=AB(3);
        vx=-uy;
        vy=ux;
        u_bg = AB(1); 
    v_bg = AB(4);
        
    %divergence and normal strain = 0 
    elseif (isDeltaZero == 1 && isZetaZero == 0 && isSigma_nZero == 1 && isSigma_sZero == 0)  
        
        R(1:n/2,2) = 0;
        R(n/2+1:n,6) = 0;
        
        R = [R(:,1) R(:,3:5)]; % remove empty columns
        
        %estimate variables
        AB = ((R'*R)\R')*UV;
        ux=0;
        uy=AB(2); % no entry for ux, so uy becomes the 2nd entry
        vx=AB(4); % no entry for ux, so vx becomes the 4th entry
        vy=0;
        u_bg = AB(1); 
    v_bg = AB(3);
        
    %vorticity and sheer strain = 0 
    elseif (isDeltaZero == 0 && isZetaZero == 1 && isSigma_nZero == 0 && isSigma_sZero == 1)  
        
        R(1:n/2,3) = 0;
        R(n/2+1:n,5) = 0;
        
        R = [R(:,1:2) R(:,4) R(:,6)]; % remove empty columns
        
        %estimate variables
        AB = ((R'*R)\R')*UV;
        ux=AB(2);
        uy=0;
        vx=0; 
        vy=AB(4); % no entry for uy or vx, so vy becomes the 4th entry  
        u_bg = AB(1); 
    v_bg = AB(3);
        
    %divergence and shear strain = 0 
    elseif (isDeltaZero == 1 && isZetaZero == 0 && isSigma_nZero == 0 && isSigma_sZero == 1)  
        
        R(n/2+1:n,6) = 0;
        R(n/2+1:n,2) = -y(:);
        R(n/2+1:n,5) = 0;
        R(n/2+1:n,3) = -x(:);
        
        R = R(:,1:4); % remove empty columns
        
        %estimate variables
        AB = ((R'*R)\R')*UV;
        ux=AB(2);
        uy=AB(3);
        vx=-uy; 
        vy=-ux;  
        u_bg = AB(1); 
    v_bg = AB(4);
        
     %vorticity and normal strain = 0 
    elseif (isDeltaZero == 0 && isZetaZero == 1 && isSigma_nZero == 1 && isSigma_sZero == 0)  
        
        R(n/2+1:n,6) = 0;
        R(n/2+1:n,2) = y(:);
        R(n/2+1:n,5) = 0;
        R(n/2+1:n,3) = x(:);
        
        R = R(:,1:4); % remove empty columns
        
        %estimate variables
        AB = ((R'*R)\R')*UV;
        ux=AB(2);
        uy=AB(3);
        vx=uy; 
        vy=ux;    
        u_bg = AB(1); 
    v_bg = AB(4);
        
        
     %%%%%%%% three variables are zero %%%%%%%% %   
        
     %divergence, vorticity and normal strain = 0 
    elseif (isDeltaZero == 1 && isZetaZero == 1 && isSigma_nZero == 1 && isSigma_sZero == 0)  
        
        R(1:n/2,2) = 0;
        R(n/2+1:n,6) = 0;
        R(n/2+1:n,5) = 0;
        R(n/2+1:n,3) = x(:);
        
        R = [R(:,1) R(:,3:4)]; % remove empty columns
        
        %estimate variables
        AB = ((R'*R)\R')*UV;
        ux=0;
        uy=AB(2); % no entry for ux, so uy becomes the 2nd entry
        vx=uy; 
        vy=0;  
        u_bg = AB(1); 
    v_bg = AB(3);
        
    %divergence, vorticity and shear strain = 0 
    elseif (isDeltaZero == 1 && isZetaZero == 1 && isSigma_nZero == 0 && isSigma_sZero == 1)  
        
        R(1:n/2,3) = 0;
        R(n/2+1:n,5) = 0;
        R(n/2+1:n,6) = 0;
        R(n/2+1:n,2) = -y(:);
        
        R = [R(:,1:2) R(:,4)]; % remove empty columns
        
        %estimate variables
        AB = ((R'*R)\R')*UV;
        ux=AB(2);
        uy=0;
        vx=0; 
        vy=-ux; 
        u_bg = AB(1); 
    v_bg = AB(3);
       
     %divergence, normal and shear strain = 0 
    elseif (isDeltaZero == 1 && isZetaZero == 0 && isSigma_nZero == 1 && isSigma_sZero == 1)  
        
        R(1:n/2,2) = 0;
        R(n/2+1:n,6) = 0;
        R(n/2+1:n,5) = 0;
        R(n/2+1:n,3) = -x(:);
        
        R = [R(:,1) R(:,3:4)]; % remove empty columns
        
        %estimate variables
        AB = ((R'*R)\R')*UV;
        ux=0;
        uy=AB(2); % no entry for ux, so uy becomes the 2nd entry
        vx=-uy; 
        vy=0;    
        u_bg = AB(1); 
    v_bg = AB(3);
        
    %vorticity, normal and shear strain = 0 
    elseif (isDeltaZero == 0 && isZetaZero == 1 && isSigma_nZero == 1 && isSigma_sZero == 1)  
        
        R(1:n/2,3) = 0;
        R(n/2+1:n,5) = 0;
        R(n/2+1:n,6) = 0;
        R(n/2+1:n,2) = y(:);
        
        R = [R(:,1:2) R(:,4)]; % remove empty columns
        
        %estimate variables
        AB = ((R'*R)\R')*UV;
        ux=AB(2);
        uy=0;
        vx=0; 
        vy=ux;  
        u_bg = AB(1); 
    v_bg = AB(3);
        
    elseif (isDeltaZero == 1 && isZetaZero ==1 && isSigma_nZero == 1 && isSigma_sZero == 1) % if all variables are zero then the solution is known and estimation isn't required.
    
        R = [R(:,1) R(:,4)];
        
        AB = ((R'*R)\R')*UV;
        ux=0;
        uy=0;
        vx=0;
        vy=0;
        u_bg = AB(1); 
    v_bg = AB(2);
        
    else
        
        msg = 'Error occurred. Unknown error.';
        error(msg)
        
    end
        
    div=ux+vy;
    vort=vx-uy;
    nstrain=ux-vy;
    sstrain=vx+uy;
    
    %Diffusivity estimation - calculate mean squared distances - must
    %calculate diffusivity outside function using moving time window -
    %can't get slope from single data point!
    
    EF = UV-R*((R'*R)\R')*UV;

    for i=1:size(u,2)
        E1(:,i) = EF(i*size(u,1)-size(u,1)+1:i*size(u,1));
        F1(:,i) = EF(length(EF)/2+i*size(u,1)-size(u,1)+1:length(EF)/2+i*size(u,1));
    end

    E=dt*cumsum(E1)';
    F=dt*cumsum(F1)';
% 
%     [M N] = size(E);
% 
%     msd_x = mean(((E'-E(:,1)' )).^2,2);
%     msd_y = mean(((F'-F(:,1)')).^2,2);
% 
%     x=1:1:N;
%     [p_x, s_x] = polyfit(x',msd_x,1);
%     f_x = polyval(p_x,1:N);
% 
%     [p_y, s_y] = polyfit(x',msd_y,1);
%     f_y = polyval(p_y,1:N);
% 
%     k_x=(f_x(N)-f_x(1))/(2*N*dt);
%     k_y=(f_y(N)-f_y(1))/(2*N*dt);
%     
%     k_x(2) = mean(dt/(2*N).*sum(E1).^2);
%     k_y(2) = mean(dt/(2*N).*sum(F1).^2);
    
    E = E1;
    F = F1;
    
    

end
