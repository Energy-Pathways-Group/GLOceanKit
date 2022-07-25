classdef TransformOperation < handle
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    properties
        name % ?
        outputVariables % Make private?
        nVarOut
        detailedDescription
        f = []
    end

    methods
        function self = TransformOperation(name,outputVariables,f)
            arguments
                name char {mustBeNonempty}
                outputVariables StateVariable {mustBeNonempty}
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

        function varargout = Compute(self,wvt,varargin)
            varargout = cell(1,self.nVarOut);
            if nargin > 0
                [varargout{:}] = self.f(wvt,varargin{:});
            else
                [varargout{:}] = self.f(wvt);
            end
        end
    end
end