function summarizeDegreesOfFreedom(self)
fprintf('----------Spatial domain----------\n');
fprintf('The spatial domain has a grid of (Nx, Ny)=(%d, %d).\n',self.Nx,self.Ny);

fprintf('The variables (u,v) each have (Nx-1)*(Ny-1)=%d degrees-of-freedom after removing the unresolved Nyquist mode.\n',(self.Nx-1)*(self.Ny-1));
fprintf('The variable eta has (Nx-1)*(Ny-1)-1=%d degrees-of-freedom after losing one degrees-of-freedom due to lack of mean eta\n',(self.Nx-1)*(self.Ny-1)-1)
fprintf('In total, this system has %d degrees-of-freedom.\n',3*(self.Nx-1)*(self.Ny-1) - 1);

fprintf('\n----------Spectral domain----------\n');
fprintf('The primary flow components have the following degrees-of-freedom:\n')
totalDOF = 0;
for name = self.flowComponentNames
    solutionGroup = self.flowComponentWithName(name);
    totalDOF = totalDOF + solutionGroup.degreesOfFreedomPerMode*solutionGroup.nModes;
    fprintf('\t%s: %d unique solutions, each with %d degrees-of-freedom.\n',solutionGroup.name,solutionGroup.nModes,solutionGroup.degreesOfFreedomPerMode);
end

fprintf('This results in a total of %d active degrees-of-freedom.\n',totalDOF);

if self.shouldAntialias == 1
    discardedModes = WVGeometryDoublyPeriodic.maskForAliasedModes(self.k_dft,self.l_dft);
    discardedModes = discardedModes & ~WVGeometryDoublyPeriodic.maskForNyquistModes(self.Nx,self.Ny);
    dof = WVGeometryDoublyPeriodic.degreesOfFreedomForRealMatrix(self.Nx,self.Ny,conjugateDimension=self.conjugateDimension);

    discardedDOFUV = sum(discardedModes(:).*dof(:));
    discardedDOFEta = sum(discardedModes(:).*dof(:));



    % discardedDOFUV = sum(discardedModes(:).*dof(:))*(self.Nj);
    % discardedDOFEta = sum(discardedModes(:).*dof(:))*(self.Nj-2);

    % discardedDOFVertical = 2*(discardedDOFUV_Nz-discardedDOFUV) + discardedDOFEta_Nz-discardedDOFEta;

    discardedDOF = 2*discardedDOFUV + discardedDOFEta;
    fprintf('There are %d modes discarded to prevent quadradic aliasing.\n',discardedDOF);
    fprintf('Active (%d) + aliased (%d) modes = %d modes\n',totalDOF,discardedDOF,discardedDOF+totalDOF);
end

% fprintf('The extra degree-of-freedom is because there is an additional constraint on the MDA modes imposed by the requirement that int N^2 eta dV=0.\n');
end