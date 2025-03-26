classdef WVFastTransformDoublyPeriodicMatlab < WVFastTransformDoublyPeriodic
    properties
        wvg

        % memory buffer to hold the dft matrices
        %
        % - Topic: Domain attributes — WV grid
        complexBuffer
    end

    methods
        function self = WVFastTransformDoublyPeriodicMatlab(wvg,Nz)
            self.wvg = wvg;
            self.complexBuffer = complex(zeros([wvg.Nx wvg.Ny Nz]));
        end
    end

    methods
        u_bar = transformFromSpatialDomainWithFourier(self,u)
        u = transformToSpatialDomainWithFourier(self,u_bar)
        du = diffX(wvg,u,options)
        du = diffY(wvg,u,options)
    end
end