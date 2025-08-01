classdef ClassB < handle
    properties (GetAccess=private, SetAccess=private)
        commonProperty
    end

    methods
        function self = ClassB
            self.commonProperty = 1;
        end
    end
end