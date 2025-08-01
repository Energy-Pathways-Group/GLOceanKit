function summarizeModeEnergy(self,options)
    % List the most energetic modes
    %
    % At the moment the +/- waves are simply added together for each mode.
    % It would be better if they were separate.
    %
    % - Topic: Energetics
    % - Declaration: summarizeModeEnergy(options)
    % - Parameter n: (optional) number of modes to list
    arguments
        self WVTransform {mustBeNonempty}
    end
    arguments
        options.n (1,1) double = 10
    end
    n = options.n;
    
    totalEnergy = self.totalEnergy;

    for name = self.flowComponentNames
        flowComponent = self.flowComponentWithName(name{1});
        flowEnergy = 0;
        if flowComponent.hasWaveComponent
            flowEnergy = flowEnergy + self.Apm_TE_factor(:).*( flowComponent.maskAp(:).*abs(self.Ap(:)).^2 + flowComponent.maskAm(:).*abs(self.Am(:)).^2 );
        end
        if flowComponent.hasPVComponent
            flowEnergy = flowEnergy +  self.A0_TE_factor(:).*( flowComponent.maskA0(:).*abs(self.A0(:)).^2);
        end
        [sortedFlowEnergy,indices] = sort(flowEnergy,'descend');

        Mode = cell(n,1);
        ConstituentEnergyPct = cell(n,1);
        OverallEnergyPct = cell(n,1);
        for iMode=1:n
            [kMode,lMode,jMode] = self.modeNumberFromIndex(indices(iMode));
            Mode{iMode} = sprintf('(%d,%d,%d)',kMode,lMode,jMode);
            ConstituentEnergyPct{iMode} = sprintf('%.3f',(sortedFlowEnergy(iMode)/sum(flowEnergy(:)))*100);
            OverallEnergyPct{iMode} = sprintf('%.3f',(sortedFlowEnergy(iMode)/totalEnergy)*100);            
        end
        Mode = string(Mode);
        ConstituentEnergyPct = string(ConstituentEnergyPct);
        OverallEnergyPct = string(OverallEnergyPct);
        if strcmp(name{1},'wave')
            Frequency = cell(n,1);
            for iMode=1:n
                Frequency{iMode} = sprintf('%.2f f',self.Omega(indices(iMode))/self.f);
            end
            Frequency = string(Frequency);
            T = table(Mode,Frequency,ConstituentEnergyPct,OverallEnergyPct);
        else
            T = table(Mode,ConstituentEnergyPct,OverallEnergyPct);
        end

        fprintf('\nThe <strong>%s</strong> flow constituent contains %.3f pct of total energy\n', self.flowComponentWithName(name{1}).name,sum(flowEnergy(:))*100/totalEnergy)
        disp(T);
    end
end