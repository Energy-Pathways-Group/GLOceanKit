function summarizeDegreesOfFreedom(self)
fprintf('Not implemented for this particular transform.\n');

return
fprintf('----------Spatial domain----------\n');
fprintf('The spatial domain has a grid of (Nx, Ny, Nz)=(%d, %d, %d).\n',self.Nx,self.Ny,self.Nz);
if self.isBarotropic == 1
    fprintf('The variables (u,v) each have (Nx-1)*(Ny-1)=%d degrees-of-freedom after removing the unresolved Nyquist mode.\n',(self.Nx-1)*(self.Ny-1));
    fprintf('The variable eta has (Nx-1)*(Ny-1)-1=%d degrees-of-freedom after losing one degrees-of-freedom due to lack of mean eta\n',(self.Nx-1)*(self.Ny-1)-1)
    fprintf('In total, this system has %d degrees-of-freedom.\n',3*(self.Nx-1)*(self.Ny-1) - 1);
else
    fprintf('The variables (u,v) each have (Nx-1)*(Ny-1)*(Nz-1)=%d degrees-of-freedom after removing the unresolved Nyquist mode.\n',(self.Nx-1)*(self.Ny-1)*(self.Nz-1));
    fprintf('The variable eta has (Nx-1)*(Ny-1)*(Nz-3)=%d degrees-of-freedom after losing two additional degrees-of-freedom due to vanishing boundaries\n',(self.Nx-1)*(self.Ny-1)*(self.Nz-3))
    fprintf('In total, this system has %d degrees-of-freedom.\n',2*(self.Nx-1)*(self.Ny-1)*(self.Nz-1) + (self.Nx-1)*(self.Ny-1)*(self.Nz-3));
end

fprintf('\n----------Spectral domain----------\n');
fprintf('The four major solutions groups have the following degrees-of-freedom:\n')
totalDOF = 0;

solutionGroup = WVGeostrophicComponent(self);
totalDOF = totalDOF + 2*solutionGroup.nModes;
fprintf('\tGeostrophic: %d unique solutions, each with 2 degrees-of-freedom.\n',solutionGroup.nModes);

solutionGroup = WVInternalGravityWaveComponent(self);
totalDOF = totalDOF + 2*solutionGroup.nModes;
fprintf('\tInternal gravity wave: %d unique solutions, each with 2 degrees-of-freedom.\n',solutionGroup.nModes);

solutionGroup = WVInertialOscillationComponent(self);
totalDOF = totalDOF + 2*solutionGroup.nModes;
fprintf('\tInertial oscillations: %d unique solutions, each with 2 degrees-of-freedom.\n',solutionGroup.nModes);

solutionGroup = WVMeanDensityAnomalyComponent(self);
totalDOF = totalDOF + solutionGroup.nModes;
fprintf('\tMean density anomaly: %d unique solutions, each with 1 degree-of-freedom.\n',solutionGroup.nModes);

fprintf('This results in a total of %d active degrees-of-freedom.\n',totalDOF);

if self.shouldAntialias == 1
    discardedModes = WVGeometryDoublyPeriodic.maskForAliasedModes(self.horizontalModes.k_dft,self.horizontalModes.l_dft);
    discardedModes = discardedModes & ~WVGeometryDoublyPeriodic.maskForNyquistModes(self.Nx,self.Ny);
    dof = WVGeometryDoublyPeriodic.degreesOfFreedomForRealMatrix(self.Nx,self.Ny,conjugateDimension=self.conjugateDimension);
    if self.isBarotropic == 1
        discardedDOFUV = sum(discardedModes(:).*dof(:));
        discardedDOFEta = sum(discardedModes(:).*dof(:));
    else
        discardedDOFUV = sum(discardedModes(:).*dof(:))*(self.Nz-1);
        discardedDOFEta = sum(discardedModes(:).*dof(:))*(self.Nz-3);
    end



    % discardedDOFUV = sum(discardedModes(:).*dof(:))*(self.Nj);
    % discardedDOFEta = sum(discardedModes(:).*dof(:))*(self.Nj-2);

    % discardedDOFVertical = 2*(discardedDOFUV_Nz-discardedDOFUV) + discardedDOFEta_Nz-discardedDOFEta;

    discardedDOF = 2*discardedDOFUV + discardedDOFEta;
    fprintf('There are %d modes discarded to prevent quadradic aliasing.\n',discardedDOF);
    fprintf('Active (%d) + aliased (%d) modes = %d modes\n',totalDOF,discardedDOF,discardedDOF+totalDOF);
end

fprintf('The extra degree-of-freedom is because there is an additional constraint on the MDA modes imposed by the requirement that int N^2 eta dV=0.\n');
end