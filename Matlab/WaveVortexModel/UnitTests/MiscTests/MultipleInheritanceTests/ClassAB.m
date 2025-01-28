classdef ClassAB < ClassA & ClassB
    properties (GetAccess=public, SetAccess=public)
        commonProperty
    end

    methods
        function self = ClassAB
            self.commonProperty = 12;
        end

        function c = foo(self,a,b)
            c = a+b;
        end
    end
end