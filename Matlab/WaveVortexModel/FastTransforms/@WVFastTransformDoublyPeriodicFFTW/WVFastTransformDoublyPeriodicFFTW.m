classdef WVFastTransformDoublyPeriodicFFTW < WVFastTransformDoublyPeriodic
    properties
        wvg

        complexBuffer

        % k wavenumber dimension on the half-complex grid
        %
        % - Topic: Domain attributes — DFT grid
        k_hc

        % l wavenumber dimension on the half-complex grid
        %
        % - Topic: Domain attributes — DFT grid
        l_hc

        % discrete Fourier transform object
        %
        % - Topic: Domain attributes — Spatial grid
        dftXY
        dftX
        dftY
        dftXYComplexBuffer
        dftXComplexBuffer
        dftYComplexBuffer
        dftRealBuffer
        Nz

        dx
        dy
    end

    methods
        function self = WVFastTransformDoublyPeriodicFFTW(wvg,Nz,options)
            arguments
                wvg 
                Nz 
                options.nCores = 12
            end
            self.wvg = wvg;
            self.Nz = Nz;
            self.complexBuffer = complex(zeros([wvg.Nx wvg.Ny Nz]));

            self.dftXY = RealToComplexTransform([wvg.Nx wvg.Ny Nz],dims=[1 2],nCores=options.nCores,planner="measure");
            self.dftX = RealToComplexTransform([wvg.Nx wvg.Ny Nz],dims=1,nCores=options.nCores,planner="measure");
            self.dftY = RealToComplexTransform([wvg.Nx wvg.Ny Nz],dims=2,nCores=options.nCores,planner="measure");

            self.dftXYComplexBuffer = complex(zeros(self.dftXY.complexSize),zeros(self.dftXY.complexSize));
            self.dftXComplexBuffer = complex(zeros(self.dftX.complexSize),zeros(self.dftX.complexSize));
            self.dftYComplexBuffer = complex(zeros(self.dftY.complexSize),zeros(self.dftY.complexSize));
            self.dftRealBuffer = double(zeros(self.dftX.realSize));

            self.k_hc = wvg.k_dft(1:(wvg.Nx/2+1));
            self.l_hc = wvg.l_dft(1:(wvg.Ny/2+1));

            self.dx = sqrt(-1)*self.k_hc;
            self.dy = sqrt(-1)*shiftdim(self.l_hc,-1);
        end
    end

    methods
        u_bar = transformFromSpatialDomainWithFourier(self,u)
        u = transformToSpatialDomainWithFourier(self,u_bar)
        du = diffX(wvg,u,n)
        du = diffY(wvg,u,n)

        du = diffXIntoArray(self,u,du,options)
    end
end