function summarizeModeEnergy(self,options)
    % List the most energetic modes
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

    A = self.Ap;
    A(1,1,:) = 0;
    C = self.Apm_TE_factor;
    energy = C(:).* (A(:).*conj(A(:)));
    [energy,indices] = sort(energy,'descend');
    totalConstituentEnergy = self.internalWaveEnergyPlus;

    fprintf('\n(k,l,j)\t\t|Ap pct\t|total pct\n');
    fprintf('----------------------------------\n');
    for iMode=1:n
        [k,l,j] = ind2sub(size(A),indices(iMode));
        fprintf('(%d,%d,%d)\t|%.3f\t|%.3f\n',k-1,l-1,j-1,(energy(iMode)/totalConstituentEnergy)*100,(energy(iMode)/totalEnergy)*100)
    end

    A = self.Am;
    A(1,1,:) = 0;
    C = self.Apm_TE_factor;
    energy = C(:).* (A(:).*conj(A(:)));
    [energy,indices] = sort(energy,'descend');
    totalConstituentEnergy = self.internalWaveEnergyMinus;

    fprintf('\n(k,l,j)\t\t|Am pct\t|total pct\n');
    fprintf('----------------------------------\n');
    for iMode=1:n
        [k,l,j] = ind2sub(size(A),indices(iMode));
        fprintf('(%d,%d,%d)\t|%.3f\t|%.3f\n',k-1,l-1,j-1,(energy(iMode)/totalConstituentEnergy)*100,(energy(iMode)/totalEnergy)*100)
    end

    A = self.A0;
    A(1,1,:) = 0;
    C = self.A0_TE_factor;
    energy = C(:).* (A(:).*conj(A(:)));
    [energy,indices] = sort(energy,'descend');
    totalConstituentEnergy = self.geostrophicEnergy;

    fprintf('\n(k,l,j)\t\t|A0 pct\t|total pct\n');
    fprintf('----------------------------------\n');
    for iMode=1:n
        [k,l,j] = ind2sub(size(A),indices(iMode));
        fprintf('(%d,%d,%d)\t|%.3f\t|%.3f\n',k-1,l-1,j-1,(energy(iMode)/totalConstituentEnergy)*100,(energy(iMode)/totalEnergy)*100)
    end

    App = self.Ap;
    Amm = self.Am;
    C = self.Apm_TE_factor;
    energy = squeeze(C(1,1,:).* (abs(App(1,1,:)).^2 + abs(Amm(1,1,:)).^2)) ;
    [energy,indices] = sort(energy,'descend');
    totalConstituentEnergy = self.inertialEnergy;

    fprintf('\n(k,l,j)\t\t|IO pct\t|total pct\n');
    fprintf('----------------------------------\n');
    for iMode=1:n
        fprintf('(1,1,%d)\t|%.3f\t|%.3f\n',indices(iMode)-1,(energy(iMode)/totalConstituentEnergy)*100,(energy(iMode)/totalEnergy)*100)
    end
end