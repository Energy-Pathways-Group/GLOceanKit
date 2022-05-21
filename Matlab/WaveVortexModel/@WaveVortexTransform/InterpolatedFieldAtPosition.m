function varargout = InterpolatedFieldAtPosition(self,x,y,z,method,varargin)
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
    dz = self.z(2)-self.z(1);
    z_index = self.Nz+floor(z/dz)-1;

    if length(unique(diff(self.z)))>1
       error('This optimized interpolation does not yet work for uneven z') 
    end

    % Identify the particles along the interpolation boundary
    if strcmp(method,'spline')
        S = 3+1; % cubic spline, plus buffer
    elseif strcmp(method,'linear')
        S = 1+1;
    end

%             % can't handle particles at the boundary
%             xrange = mod(((min(x_index)-S):(max(x_index)+S)),self.Nx)+1;
%             yrange = mod(((min(y_index)-S):(max(y_index)+S)),self.Ny)+1;
    xrange = ((min(x_index)-S):(max(x_index)+S))+1;
    yrange = ((min(y_index)-S):(max(y_index)+S))+1;
    zrange = (max((min(z_index)-S),0):min((max(z_index)+S),self.Nz-1))+1;
    
    [X,Y,Z] = ndgrid(self.x,self.y,self.z);
    varargout = cell(1,nargout);
    for i = 1:nargout
        U = varargin{i}; % gridded field
        u = interpn(X(xrange,yrange,zrange),Y(xrange,yrange,zrange),Z(xrange,yrange,zrange),U(xrange,yrange,zrange),x_tilde,y_tilde,z,method);
        varargout{i} = u;
    end
end