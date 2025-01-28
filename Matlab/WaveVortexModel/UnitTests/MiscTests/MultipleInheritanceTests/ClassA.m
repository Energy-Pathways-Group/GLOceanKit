classdef ClassA < handle
    properties (GetAccess=private, SetAccess=private)
        commonProperty
    end

    methods
        function self = ClassA
            self.commonProperty = 2;
        end

        function c = foo(self,a)
            c = a+a;
        end
    end
end