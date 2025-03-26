classdef WVFastTransformDoublyPeriodicFFTW < WVFastTransformDoublyPeriodic
    properties
        wvg

        % memory buffer to hold the dft matrices
        %
        % - Topic: Domain attributes — WV grid
        complexBuffer

        fftw_complex_cache
        fftw_real_cache


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
        dft
    end

    methods
        function self = WVFastTransformDoublyPeriodicFFTW(wvg)
            self.wvg = wvg;
            self.complexBuffer = complex(zeros([wvg.Nx wvg.Ny wvg.Nz]));
            self.dft = RealToComplexTransform([Nxy(1) Nxy(2) options.Nz],dims=[1 2],nCores=8,planner="measure");
            self.fftw_complex_cache = complex(zeros(self.dft.complexSize));
            self.fftw_real_cache = double(zeros(self.dft.complexSize));
            self.k_hc = wvg.k_dft(1:(wvg.Nx/2+1));
            self.l_hc = wvg.l_dft(1:(wvg.Ny/2+1));
        end
    end

    methods
        u_bar = transformFromSpatialDomainWithFourier(self,u)
        u = transformToSpatialDomainWithFourier(self,u_bar)
        du = diffX(wvg,u,n)
        du = diffY(wvg,u,n)
    end
end