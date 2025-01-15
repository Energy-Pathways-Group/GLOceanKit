classdef ClassA < handle
    properties (GetAccess=private, SetAccess=private)
        commonProperty
    end

    methods
        function self = ClassA
            self.commonProperty = 2;
        end
    end
end