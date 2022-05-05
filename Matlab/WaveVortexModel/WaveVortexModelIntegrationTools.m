classdef WaveVortexModelIntegrationTools < handle
    %WaveVortexModelIntegrationTools Tools for integrating (time-stepping)
    %the WaveVortexModel.
    %
    %   inttool = WaveVortexModelIntegrationTools(wvm,newModelOutput,outputInterval) creates a new
    %   integration tool for the model. newModelOutput is optional (you
    %   don't have to output to a file).
    %
    %   inttool = WaveVortexModelIntegrationTools(existingModelOutput,newModelOutput,outputInterval)
    %   opens existing NetCDF output from the WaveVortexModel and uses that
    %   for a restart. newModelOutput is optional (you don't have to output
    %   to a file).
    %
    %   'shouldOverwriteExisting', default 0
    %   'shouldDoubleResolution', 0 or 1
    %   'restartIndex', index in existingModelOutput to use as restart.
    properties
        wvm             % WaveVortexModel
        t=0             % current model time (in seconds)
    
        outputFile
        netcdfTool      % WaveVortexModelNetCDFTools instance---empty indicates no file output
        outputInterval  % model output interval (seconds)
        t0=0            % initial output time
        iTime=1         % current index of outputTimes written to file.

        integrator      % Array integrator
        
        isDynamicsLinear = 0
        xFloat, yFloat, zFloat
        xDrifter, yDrifter, zDrifter
        tracers
    end

    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = WaveVortexModelIntegrationTools(varargin)
            newModelOutput = [];
            if nargin > 1
                newModelOutput = varargin{2};
            end
            outputInterval = [];
            if nargin > 2
                outputInterval = varargin{3};
            end
            existingModelOutput = [];

            extraargs = varargin(4:end);
            if mod(length(extraargs),2) ~= 0
                error('Arguments must be given as name/value pairs.');
            end
            shouldOverwriteExisting = 0;

            if isa(varargin{1},'WaveVortexModel')
                waveVortexModel = varargin{1};
                for k = 1:2:length(extraargs)
                    if strcmp(extraargs{k}, 'shouldOverwriteExisting')
                        shouldOverwriteExisting = extraargs{k+1};
                    else
                        error('Unknown argument, %s', extraargs{k});
                    end
                end
            elseif isa(varargin{1},'char' )
                existingModelOutput = varargin{1};

                shouldDoubleResolution = 0;
                restartIndex = Inf;
                for k = 1:2:length(extraargs)
                    if strcmp(extraargs{k}, 'shouldDoubleResolution')
                        shouldDoubleResolution = extraargs{k+1};
                    elseif strcmp(extraargs{k}, 'restartIndex')
                        restartIndex = extraargs{k+1};
                    elseif strcmp(extraargs{k}, 'shouldOverwriteExisting')
                        shouldOverwriteExisting = extraargs{k+1};
                    else
                        error('Unknown argument, %s', extraargs{k});
                    end
                end

                nctool = WaveVortexModelNetCDFTools(existingModelOutput,'timeIndex',restartIndex);

                self.t = nctool.t;
                self.t0 = self.t;
                if (shouldDoubleResolution == 0)
                    waveVortexModel = nctool.wvm;
                else
                    waveVortexModel = nctool.wvm.waveVortexModelWithResolution(2*[nctool.wvm.Nx,nctool.wvm.Ny,nctool.wvm.nModes]);
                end
            end

            self.wvm = waveVortexModel;
            self.outputFile = newModelOutput;
            if ~isempty(self.outputFile)
                % user wants us to output to a netcdf file
                
                % Figure out the output interval
                if ~isempty(outputInterval)
                    % if they set an output interval, obviously use that
                    self.outputInterval = outputInterval;
                elseif ~isempty(existingModelOutput)
                    % if there's existing model output, use that output
                    % interval
                    time = ncread(existingModelOutput,'t');
                    if length(time)>1
                        self.outputInterval = time(2)-time(1);
                    else
                        self.outputInterval = self.wvm.inertialPeriod/10;
                    end
                else
                    % otherwise hit our default
                    self.outputInterval = self.wvm.inertialPeriod/10;
                end

                if ~isempty(existingModelOutput) && strcmp(self.outputFile,existingModelOutput) == 1
                    % okay, they want us to append to the existing file.
                    % We just need to set t0 and iTime accordingly
                    error('Cannot append to existing file!');
                    % Can't do this yet because the NetCDFTools don't fetch
                    % all the varIDs, etc.
                    self.netcdfTool = nctool;
                    time = ncread(nctool.netcdfFile,'t');
                    self.t0 = time(1);
                    self.iTime = length(time);
                else
                    [filepath,name,~] = fileparts(self.outputFile);
                    matFilePath = sprintf('%s/%s.mat',filepath,name);
                    if isfile(self.outputFile) || isfile(matFilePath)
                        if shouldOverwriteExisting == 1
                            if isfile(self.outputFile)
                                delete(self.outputFile);
                            end
                            if isfile(matFilePath)
                                delete(matFilePath);
                            end
                        else
                            error('File already exists!');
                        end
                    end
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Integration loop
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = IntegrateToTime(self,time,cfl)
            if nargin < 3 || isempty(cfl)
                cfl=0.5;
            end
            self.SetupIntegrator(self.t,cfl);
            
            % total dT time steps to meet or exceed the requested time.
            nIncrements = ceil((time-self.t)/self.integrator.stepSize);
            incrementsWrittenToFile = 0;

            if ~isempty(self.outputFile)
                self.SetupNetCDFToolsForOutputFile();
                stepsPerOutput = round(self.outputInterval/self.integrator.stepSize);
                firstIncrement = round((self.t0-self.t)/self.integrator.stepSize);
                outputIncrements = firstIncrement:stepsPerOutput:nIncrements;
                outputTimes = self.t + outputIncrements*self.integrator.stepSize;

                % Save the initial conditions
                if outputIncrements(self.iTime)==0 
                    self.netcdfTool.WriteTimeAtIndex(self.iTime,self.t);
                    self.netcdfTool.WriteAmplitudeCoefficientsAtIndex(self.iTime);
                    self.netcdfTool.WriteEnergeticsAtIndex(self.iTime);
                    self.netcdfTool.WriteEnergeticsKJAtIndex(self.iTime);
                    self.iTime = self.iTime + 1;
                    incrementsWrittenToFile = incrementsWrittenToFile + 1;
                end
            end
            
            for i=1:nIncrements
                if i == 1
                    fprintf('Starting numerical simulation on %s.\n', datestr(datetime('now')));
                    inertialPeriod = (2*pi/(2 * 7.2921E-5 * sin( self.wvm.latitude*pi/180 )));
                    fprintf('\tStarting at model time t=%.2f inertial periods and integrating to t=%.2f inertial periods with %d RK4 time steps.\n',self.t/inertialPeriod,time/inertialPeriod,nIncrements);
                    if ~isempty(self.outputFile)
                        fprintf('\tWriting %d of those time steps to file. Will write to output file starting at index %d.\n',sum(outputIncrements>=0),self.iTime-incrementsWrittenToFile);
                    end
                elseif i == 2
                    startTime = datetime('now');
                else
                    timePerStep = (datetime('now')-startTime)/(i-2);
                    % We want to inform the user about every 30 seconds
                    stepsPerInform = ceil(30/seconds(timePerStep));
                    if (i==3 || mod(i,stepsPerInform) == 0)
                        timeRemaining = (nIncrements-i+1)*timePerStep;
                        fprintf('\tmodel time t=%.2f inertial periods, RK4 time step %d of %d. Estimated finish time %s (%s from now)\n', self.t/inertialPeriod, i, nIncrements, datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
                        self.wvm.summarizeEnergyContent();
                    end
                end
                
                % TODO: update this handle variable time-stepped inputs
                self.integrator.IncrementForward();
                self.wvm.Ap = self.integrator.currentY{1};
                self.wvm.Am = self.integrator.currentY{2};
                self.wvm.A0 = self.integrator.currentY{3};
                self.t = self.integrator.currentTime;

                if ~isempty(self.outputFile) && self.iTime <= length(outputIncrements)
                    if outputIncrements(self.iTime) == i   
                        self.netcdfTool.WriteTimeAtIndex(self.iTime,outputTimes(self.iTime));
                        self.netcdfTool.WriteAmplitudeCoefficientsAtIndex(self.iTime);
                        self.netcdfTool.WriteEnergeticsAtIndex(self.iTime);
                        self.netcdfTool.WriteEnergeticsKJAtIndex(self.iTime);
                        self.netcdfTool.sync();
                        self.iTime = self.iTime + 1;
                        incrementsWrittenToFile = incrementsWrittenToFile + 1;
                    end
                end
            end

            if ~isempty(self.outputFile)
                fprintf('Ending simulation. Wrote %d time points to file\n',incrementsWrittenToFile);
                self.netcdfTool.close();
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Internal
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = SetupNetCDFToolsForOutputFile(self)
            if isempty(self.outputFile)
                % user didn't request output, so move on
                return;
            end

            if isempty(self.netcdfTool)
                self.netcdfTool = WaveVortexModelNetCDFTools(self.wvm,self.outputFile);
                self.netcdfTool.CreateAmplitudeCoefficientVariables();
                self.netcdfTool.CreateEnergeticsVariables();
                self.netcdfTool.CreateEnergeticsKJVariables();
            else
                if isempty(self.netcdfTool.ncid)
                    self.netcdfTool.open();
                end
            end

            
        end

        function self = SetupIntegrator(self,initialTime,cfl)
            deltaT = self.wvm.TimeStepForCFL(cfl,self.outputInterval);
            Y0 = self.InitialConditionsArray();
            if isempty(Y0{1})
                error('Nothing to do! You must have set to linear dynamics, without floats, drifters or tracers.');
            end

            self.integrator = ArrayIntegrator(@(t,y0) self.FluxAtTime(t,y0),Y0,deltaT);
            self.integrator.currentTime = initialTime;
        end
        
        function Y0 = InitialConditionsArray(self)
            Y0 = cell(1,1);
            n = 0;
            if self.isDynamicsLinear == 0
                n=n+1;Y0{n} = self.wvm.Ap;
                n=n+1;Y0{n} = self.wvm.Am;
                n=n+1;Y0{n} = self.wvm.A0;
            end

            if ~isempty(self.xFloat)
                n=n+1;Y0{n} = self.xFloat;
                n=n+1;Y0{n} = self.yFloat;
                n=n+1;Y0{n} = self.zFloat;
            end

            if ~isempty(self.xDrifter)
                n=n+1;Y0{n} = self.xDrifter;
                n=n+1;Y0{n} = self.yDrifter;
            end

            if ~isempty(self.tracers)
                for i=1:length(self.tracers)
                    n=n+1;Y0{n} = self.tracers{i};
                end
            end
        end

        function F = FluxAtTime(self,t,y0)
            F = cell(1,1);
            n = 0;
            if self.isDynamicsLinear == 0
                [Fp,Fm,F0,U,V,W] = self.wvm.NonlinearFluxAtTime(t,y0{n+1},y0{n+2},y0{n+3});
                n=n+1;F{n} = Fp;
                n=n+1;F{n} = Fm;
                n=n+1;F{n} = F0;
            else
                [U,V,W] = self.wvm.VelocityFieldAtTime(t);
            end

            if ~isempty(self.xFloat)
                [Fx,Fy,Fz] = self.InterpolatedFieldAtPosition(y0{n+1},y0{n+2},y0{n+3},'spline',U,V,W);
                n=n+1;Y0{n} = Fx;
                n=n+1;Y0{n} = Fy;
                n=n+1;Y0{n} = Fz;
            end

            if ~isempty(self.xDrifter)
                [Fx,Fy] = self.InterpolatedFieldAtPosition(y0{n+1},y0{n+2},self.zDrifter,'spline',U,V);
                n=n+1;Y0{n} = Fx;
                n=n+1;Y0{n} = Fy;
            end

            if ~isempty(self.tracers)
                for i=1:length(self.tracers)
                    phibar = self.TransformFromSpatialDomainWithF(y0{n+1});
                    [~,Phi_x,Phi_y,Phi_z] = self.TransformToSpatialDomainWithFAllDerivatives(phibar);
                    n=n+1;Y0{n} = -U.*Phi_x - V.*Phi_y - W.*Phi_z;
                end
            end
        end

    end
end