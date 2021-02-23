function [varargout] = InternalVariablesAtTimePosition(self,t,x,y,z,method,varargin)
    % Returns the gridded/internal dynamical variables at any point in
    % space or time. Because these waves are gridded, we can use
    % interpolation (e.g, spline, linear) to deduce the variable
    % values off grid, or we can use spectral interpolation to get
    % 'exact' values.
    if self.advectionSanityCheck == 0
        self.advectionSanityCheck = 1;
        if (self.z(end)-self.z(1)) ~= self.Lz
            warning('Vertical domain does not span the full depth of the ocean. This will lead to NaNs when advected particles leave the resolved domain.')
        end
    end

    varargout = cell(size(varargin));
    if strcmp(method,'exact')
        [varargout{:}] = self.InternalVariablesAtTimePositionExact(t,x,y,z,varargin{:});
    else
        [varargout{:}] = self.InternalVariableFieldsAtTime(t,varargin{:});
        [varargout{:}] = self.InterpolatedFieldAtPosition(x,y,z,method,varargout{:});
    end
end