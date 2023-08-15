function varargout = interpolatedFieldAtPosition(self,x,y,z,method,varargin)
    if nargin-5 ~= nargout
        error('You must have the same number of input variables as output variables');
    end

    % (x,y) are periodic for the gridded solution
    x_tilde = mod(x,self.Lx);
    y_tilde = mod(y,self.Ly);

    dx = self.x(2)-self.x(1);
    x_index = floor(x_tilde/dx);
    dy = self.y(2)-self.y(1);
    y_index = floor(y_tilde/dy);

    if(strcmp(method, "finufft"))
        varargout = cell(1,nargout);
        for i = 1:nargout
            U = varargin{i};
            ubar = self.transformFromSpatialDomainWithF(U);
            u = self.transformToSpatialDomainWithFNU(ubar, x_tilde, y_tilde);
            % u = self.transformToSpatialDomainWithF(ubar);
            varargout{i} = u;
        end
        return
    end

    % Identify the particles along the interpolation boundary
    if strcmp(method,'spline')
        S = 3+1; % cubic spline, plus buffer
    elseif strcmp(method,'linear')
        S = 1+1;
    end
    bpx = x_index < S-1 | x_index > self.Nx-S;
    bpy = y_index < S-1 | y_index > self.Ny-S;

    varargout = cell(1,nargout);
    for i = 1:nargout
        U = varargin{i}; % gridded field
        u = zeros(size(x)); % interpolated value
        u(~bpx & ~bpy) = interpn(self.X,self.Y,U,x_tilde(~bpx & ~bpy),y_tilde(~bpx & ~bpy),method);
        if any(bpx & bpy)
            % then do the same for particles that along the boundary.
            x_tildeS = mod(x(bpx & bpy)+S*dx,self.Lx);
            y_tildeS = mod(y(bpx & bpy)+S*dy,self.Ly);
            u(bpx & bpy) = interpn(self.X,self.Y,circshift(U,[S S 0]),x_tildeS,y_tildeS,method);
        end
        
        if any(bpx & ~bpy)
            x_tildeS = mod(x(bpx & ~bpy)+S*dx,self.Lx);
            u(bpx & ~bpy) = interpn(self.X,self.Y,circshift(U,[S 0 0]),x_tildeS,y_tilde(bpx & ~bpy),method);
        end

        if any(~bpx & bpy)
            y_tildeS = mod(y(~bpx & bpy)+S*dy,self.Ly);
            u(~bpx & bpy) = interpn(self.X,self.Y,circshift(U,[0 S 0]),x_tilde(~bpx & bpy),y_tildeS,method);
        end
        varargout{i} = u;
        if any(isnan(u(:)))
            fprintf('bad apple.\n')
        end
    end
end