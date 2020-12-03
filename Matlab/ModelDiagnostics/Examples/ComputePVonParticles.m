% Point our diagnostics class to some Winters model output
file = '/Volumes/Samsung_T5/nsf_iwv/WintersNonlinear/EarlyV2_GM_NL_forced_damped_restart';
diag = ModelDiagnosticsWintersModel(file);

% Let's load up the first saved file...
diag.setTimeIndex(1);



% grab the particle trajectories from the output
[t,x,y,z] = diag.wintersmodel.ParticleTrajectories;

% we need to check and see if the saved Eulerian files are output at
% the same time points. So grab the time-stamps from the saved files...
t_fields = zeros(diag.nT,1);
for iT = 1:diag.nT
   t_fields(iT) = diag.wintersmodel.VariableFieldsFrom3DOutputFileAtIndex(iT,'t');
end

% ...find the common set of time points
[~,validParticleIndices,validFileIndices] = intersect(t,t_fields);

% and then reduce the trajectories down to those points
t = t(validParticleIndices);
x = x(validParticleIndices,:);
y = y(validParticleIndices,:);
z = z(validParticleIndices,:);

ErtelPV = zeros(size(x));
LinearErtelPV = zeros(size(x));
QGPV = zeros(size(x));

startTime = datetime('now');
for iTimeIndex = 1:length(validFileIndices)
    if iTimeIndex>=2
        timePerStep = (datetime('now')-startTime)/(iTimeIndex-1);
        timeRemaining = (length(validFileIndices)-iTimeIndex+1)*timePerStep;
        fprintf('\ttime step %d of %d. Estimated finish time %s (%s from now)\n', iTimeIndex, length(validFileIndices), datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
    end
    
    % Load up the dynamical variables at the given time
    diag.setTimeIndex(validFileIndices(iTimeIndex));
    
    % ...and then compute some flavors of PV.
    noBackground = 1; % okay if N2=constant
    % subtracting off Lz in z to make Winters output work!!! need a better
    % fix.
    [ErtelPV(iTimeIndex,:),LinearErtelPV(iTimeIndex,:),QGPV(iTimeIndex,:)] = diag.InterpolatedFieldAtPosition(x(iTimeIndex,:),y(iTimeIndex,:),z(iTimeIndex,:)-diag.Lz,'linear',diag.ErtelPV(noBackground), diag.LinearErtelPV(noBackground), diag.QGPV);   
end

% save('ParticlePV.mat','t','x','y','z','ErtelPV','LinearErtelPV','QGPV');