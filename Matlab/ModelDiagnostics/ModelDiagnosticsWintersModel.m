classdef ModelDiagnosticsWintersModel < ModelDiagnostics
    % ModelDiagnostics
    properties (Access = public)
        wavemodel
        wintersmodel
        nT
        timeIndex
    end
    
    methods
        function self = ModelDiagnosticsWintersModel(file)
            wintersmodel_ = WintersModel(file);
            wavemodel_ = wintersmodel_.wavemodel;
            
            dims = [wavemodel_.Lx,wavemodel_.Ly,wavemodel_.Lz];
            n = [wavemodel_.Nx,wavemodel_.Ny,wavemodel_.Nz];
            latitude = wavemodel_.latitude;
            rhoBar = wavemodel_.RhoBarAtDepth(wavemodel_.z);
            N2 = wavemodel_.N2AtDepth(wavemodel_.z);
            
            self@ModelDiagnostics(dims, n, latitude, rhoBar, N2)
            
            self.wintersmodel = wintersmodel_;
            self.wavemodel = wintersmodel_;
            self.nT = self.wintersmodel.NumberOf3DOutputFiles;
        end
         
        function self = setTimeIndex(self,iT)
            [t0,u,v,w,rho_prime] = self.wintersmodel.VariableFieldsFrom3DOutputFileAtIndex(iT,'t','u','v','w','rho_prime');
            self.InitializeWithHorizontalVelocityAndDensityPerturbationFields(u,v,w,rho_prime);
            self.timeIndex = iT;
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % 
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
         
        
    end
end