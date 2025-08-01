classdef WVOperation < handle
%Perform an operation and return a variable using a WVTransform
% 
% A WVOperation allows you to add new functionality to the WVTransform by
% defining one or more variables and creating an operation for computing
% those variables.
% 
% If the operation is simple and can be computed in one line, you can
% directly instantiate the WVOperation class and pass a function handle.
% For example,
% 
% ```matlab
% outputVar = WVVariableAnnotation('zeta_z',{'x','y','z'},'1/s^2', 'vertical component of relative vorticity');
% f = @(wvt) wvt.diffX(wvt.v) - wvt.diffY(wvt.u);
% wvt.addOperation(WVOperation('zeta_z',outputVar,f));
% ```
% 
% will enable direct calls to `wvt.zeta_z` to compute the vertical
% vorticity.
% 
% Note that a `WVOperation` that computes a single variable, must have the
% same name as the variable, as specified in `WVVariableAnnotation`.
%
% More involved calculations may require subclassing WVOperation and
% overriding the `compute` method. Note that every variable returned by the
% compute operation must be described with a `WVVariableAnnotation`.
% 
% - Declaration: classdef WVOperation < handle
    properties
        % name of the operation
        %
        % If the operation only returns a single variable, it **must** be
        % given the same name as that variable.
        % - Topic: Properties
        name


        detailedDescription
    end

    properties (GetAccess=public, SetAccess=protected)
        % array of WVVariableAnnotations describing the outputs of the computation
        %
        % This array is set during initialization.
        %
        % The order of the WVVariableAnnotation must match the order that
        % the variables will be returned with the -compute method.
        % - Topic: Properties
        outputVariables

        % number of variables returned by the computation
        %
        % - Topic: Properties
        nVarOut

        % function handle to be called when computing the operation
        %
        % If you override the compute method, there no need to set or use
        % this function handle.
        % - Topic: Properties
        f = []
    end

    methods
        function self = WVOperation(name,outputVariables,f)
            % create a new WVOperation for computing a new variable
            %
            % - Topic: Initialization
            % - Declaration: operation = WVOperation(name,outputVariables,f)
            % - Parameter name: name of the operation
            % - Parameter outputVariables: ordered array of WVVariableAnnotations
            % - Parameter f: function handle that takes a WVTransform as an argument and returns variables matching outputVariable.
            % - Return operation: a new WV operation instance
            arguments
                name char {mustBeNonempty}
                outputVariables WVVariableAnnotation {mustBeNonempty}
                f function_handle
            end

            self.outputVariables = outputVariables;
            self.nVarOut = length(self.outputVariables);

            % Make each output variable aware of the operation that
            % computes it.
            for iVar=1:length(self.outputVariables)
                self.outputVariables(iVar).modelOp = self;
            end
            
            self.name = name;
            if self.nVarOut == 1
                if ~strcmp(self.name, self.outputVariables(1).name)
                    error('An operation with only one output variable must have the same name as the first output variable.')
                end
            end
            
            self.f = f;
        end

        function varargout = compute(self,wvt,varargin)
            % compute the promised variable
            %
            % The compute operation is given the current state of the ocean
            % (wvt) during each call, and it is expected to compute the
            % variables from that state. You can, of course, use other 
            % instance variables from your own custom subclass to
            % implement computations with other dependencies.
            % 
            % - Topic: Computation
            % - Declaration: varargout = compute(wvt,varargin)
            % - Parameter wvt: A WVTransform instance from which to compute the variable
            % - Returns varargout: cell array of returned variables
            varargout = cell(1,self.nVarOut);
            if ~isempty(varargin)
                [varargout{:}] = self.f(wvt,varargin{:});
            else
                [varargout{:}] = self.f(wvt);
            end
        end
    end
end