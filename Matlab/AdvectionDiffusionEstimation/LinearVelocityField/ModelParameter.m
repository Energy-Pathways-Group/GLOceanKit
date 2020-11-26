classdef ModelParameter
    enumeration
        u0v0, u1v1, strain, vorticity, divergence
    end
    methods (Static)
        function s = name(p)
            switch(p)
                case ModelParameter.u0v0
                    s = 'u_0v_0';
                case ModelParameter.u1v1
                    s = 'u_1v_1';
                case ModelParameter.strain
                    s = '\sigma';
                case ModelParameter.vorticity
                    s = '\zeta';
                case ModelParameter.divergence
                    s = '\delta';
            end
        end
        
        function s = modelName(p)
            s = '(';
            for i=1:length(p)
                s = strcat(s,ModelParameter.name(p(i)));
                if i<length(p)
                    s = strcat(s,',');
                end
            end
            s = strcat(s,')');
        end
    end
end

