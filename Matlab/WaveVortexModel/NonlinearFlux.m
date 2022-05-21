classdef NonlinearFlux < ModelOperation
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here

    methods
        function self = NonlinearFlux()
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            obj.Property1 = inputArg1 + inputArg2;
        end

        function varargout = Compute(self,wvt)
            u = wvt.u;
            v = wvt.v;
            w = wvt.w;
            eta = wvt.eta;

            uNL = u.*DiffFourier(wvt.x,u,1,1) + v.*DiffFourier(wvt.y,u,1,2) + w.*DiffCosine(wvt.z,u,1,3);
            vNL = u.*DiffFourier(wvt.x,v,1,1) + v.*DiffFourier(wvt.y,v,1,2) + w.*DiffCosine(wvt.z,v,1,3);
            nNL = u.*DiffFourier(wvt.x,eta,1,1) + v.*DiffFourier(wvt.y,eta,1,2) + w.*(DiffSine(wvt.z,eta,1,3) + eta.*wvt.dLnN2);


        end
    end
end