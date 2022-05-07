classdef WaveVortexModelIntegrationTools < handle
    %WaveVortexModelIntegrationTools Tools for integrating (time-stepping)
    %the WaveVortexModel.
    %
    %   inttool = WaveVortexModelIntegrationTools(wvm) creates a new
    %   integration tool for the model.
    %
    %   inttool = WaveVortexModelIntegrationTools(existingModelOutput,restartIndex,shouldDoubleResolution)
    %   opens existing NetCDF output from the WaveVortexModel and uses that
    %   for a restart. restartIndex is optional, defaults to Inf (last time
    %   point). shouldDoubleResolution is optional, defaults to 0.
    %
    %   'shouldOverwriteExisting', default 0
    %   'shouldDoubleResolution', 0 or 1
    %   'restartIndex', index in existingModelOutput to use as restart.

    % TODO, May 5th, 2022
    % - add multiple tracers, store ids and names in struct?
    % - remove tracer variance at aliased wavenumbers?
    % - add scalar option for floats and drifters, e.g., save density or pv
    % - linear dynamics should only save the coefficients once (actually
    %   option should be to only write initial conditions)
    % - need method of doing fancy stuff during the integration loop
    % - want to write float and drifter paths to memory
    % - Maybe a list of variables (as enums) that we want to write
    % - 
    properties
        wvm             % WaveVortexModel
        t=0             % current model time (in seconds)
    
        outputFile
        netcdfTool      % WaveVortexModelNetCDFTools instance---empty indicates no file output
        outputInterval  % model output interval (seconds)
        t0=0            % initial output time
        timeIndex=1     % current index of outputTimes written to file.
        outputIncrements
        incrementsWrittenToFile

        saveFloatsAndDriftersToMemory = 0   % set to true if no output file is specified
        xFloatT, yFloatT, zFloatT, xDrifterT, yDrifterT, zDrifterT

        integrator      % Array integrator
        startTime
        nIncrements
        
        linearDynamics = 0
        xFloat, yFloat, zFloat
        xDrifter, yDrifter, zDrifter
        tracers, tracerNames
    end

    methods

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Initialization
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function self = WaveVortexModelIntegrationTools(varargin)

            if isa(varargin{1},'WaveVortexModel')
                waveVortexModel = varargin{1};
            elseif isa(varargin{1},'char' )
                existingModelOutput = varargin{1};

                restartIndex = Inf;
                shouldDoubleResolution = 0;
                if nargin > 1
                    restartIndex = varargin{2};
                end
                if nargin > 2
                    restartIndex = varargin{3};
                end

                nctool = WaveVortexModelNetCDFTools(existingModelOutput,'timeIndex',restartIndex);

                self.t = nctool.t;
                self.t0 = self.t;
                if (shouldDoubleResolution == 0)
                    waveVortexModel = nctool.wvm;
                else
                    waveVortexModel = nctool.wvm.waveVortexModelWithResolution(2*[nctool.wvm.Nx,nctool.wvm.Ny,nctool.wvm.nModes]);
                end

                % if there's existing model output, use that output interval
                time = ncread(existingModelOutput,'t');
                if length(time)>1
                    self.outputInterval = time(2)-time(1);
                end
            end

            self.wvm = waveVortexModel;  
        end
        


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Floats and drifters and tracer!
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function SetFloatPositions(self,x,y,z)
            self.xFloat = reshape(x,1,[]);
            self.yFloat = reshape(y,1,[]);
            self.zFloat = reshape(z,1,[]);
            if (length(self.xFloat) ~= length(self.yFloat)) || (length(self.xFloat) ~= length(self.zFloat))
                error('SetFloatPositions failed! (x,y,z) must have the same length.')
            end
        end

        function SetDrifterPositions(self,x,y,z)
            self.xDrifter = reshape(x,1,[]);
            self.yDrifter = reshape(y,1,[]);
            self.zDrifter = reshape(z,1,[]);
            if (length(self.xDrifter) ~= length(self.yDrifter)) || (length(self.xDrifter) ~= length(self.zDrifter))
                error('SetDrifterPositions failed! (x,y,z) must have the same length.')
            end
        end

        function AddTracer(self,phi,name)
            if isempty(self.tracers)
                self.tracers{1} = phi;
                self.tracerNames{1} = name;
            else
                self.tracers{end+1} = phi;
                self.tracerNames{end+1} = name;
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
            self.nIncrements = ceil((time-self.t)/self.integrator.stepSize);
            
            self.OpenNetCDFFileForTimeStepping();
            
            for iIncrement=1:self.nIncrements
                self.ShowIntegrationTimeDiagnostics(iIncrement);
                
                self.integrator.IncrementForward();
                n=0;
                if self.linearDynamics == 0
                    n=n+1; self.wvm.Ap = self.integrator.currentY{n};
                    n=n+1; self.wvm.Am = self.integrator.currentY{n};
                    n=n+1; self.wvm.A0 = self.integrator.currentY{n};
                end

                if ~isempty(self.xFloat)
                    n=n+1; self.xFloat = self.integrator.currentY{n};
                    n=n+1; self.yFloat = self.integrator.currentY{n};
                    n=n+1; self.zFloat = self.integrator.currentY{n};
                end

                if ~isempty(self.xDrifter)
                    n=n+1; self.xDrifter = self.integrator.currentY{n};
                    n=n+1; self.yDrifter = self.integrator.currentY{n};
                end

                if ~isempty(self.tracers)
                    for iTracer=1:length(self.tracers)
                        n=n+1; self.tracers{iTracer} = self.integrator.currentY{n};
                    end
                end
                
                self.t = self.integrator.currentTime;

                self.WriteTimeStepToNetCDFFile(iIncrement);
            end

            self.CloseNetCDFFile();
        end

        function ShowIntegrationTimeDiagnostics(self,integratorIncrement)
            if integratorIncrement == 1
                fprintf('Starting numerical simulation on %s.\n', datestr(datetime('now')));
                fprintf('\tStarting at model time t=%.2f inertial periods and integrating to t=%.2f inertial periods with %d RK4 time steps.\n',self.t/self.wvm.inertialPeriod,time/self.wvm.inertialPeriod,self.nIncrements);
                if ~isempty(self.outputFile)
                    fprintf('\tWriting %d of those time steps to file. Will write to output file starting at index %d.\n',sum(self.outputIncrements>=0),self.timeIndex-self.incrementsWrittenToFile);
                end
            elseif integratorIncrement == 2
                self.startTime = datetime('now');
            else
                timePerStep = (datetime('now')-self.startTime)/(integratorIncrement-2);
                % We want to inform the user about every 30 seconds
                stepsPerInform = ceil(30/seconds(timePerStep));
                if (integratorIncrement==3 || mod(integratorIncrement,stepsPerInform) == 0)
                    timeRemaining = (self.nIncrements-integratorIncrement+1)*timePerStep;
                    fprintf('\tmodel time t=%.2f inertial periods, RK4 time step %d of %d. Estimated finish time %s (%s from now)\n', self.t/inertialPeriod, integratorIncrement, self.nIncrements, datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
                    self.wvm.summarizeEnergyContent();
                end
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % NetCDF Output
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function CreateNetCDFFileForModelOutput(self,modelOutputFile,outputInterval,overwriteExisting)
            %CreateNetCDFFileForModelOutput
            % modelOutputFile   file path for output
            % outputInterval    (optional) interval (in seconds) to write
            %                   to file. If restarting, will use previous time step,
            %                   otherwise will write at inertialPeriod/10.
            % overwriteExisting (optional) default to 0
            self.outputFile = modelOutputFile;

            % Figure out the output interval
            if ~isempty(outputInterval)
                % if they set an output interval, obviously use that
                self.outputInterval = outputInterval;
            elseif isempty(self.outputInterval)
                % if one hasn't been set, use the default
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
                self.timeIndex = length(time);
            else
                [filepath,name,~] = fileparts(self.outputFile);
                matFilePath = sprintf('%s/%s.mat',filepath,name);
                if isfile(self.outputFile) || isfile(matFilePath)
                    if overwriteExisting == 1
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

        function OpenNetCDFFileForTimeStepping(self)
            if isempty(self.outputFile)
                % user didn't request output, so move on
                return;
            end

            if isempty(self.netcdfTool)
                % initialize the file, and appropriate variables
                self.netcdfTool = WaveVortexModelNetCDFTools(self.wvm,self.outputFile);

                if self.linearDynamics == 0
                    self.netcdfTool.InitializeAmplitudeCoefficientStorage();
                    self.netcdfTool.InitializeEnergeticsStorage();
                    self.netcdfTool.InitializeEnergeticsKJStorage();
                end

                if ~isempty(self.xFloat)
                    self.netcdfTool.InitializeFloatStorage(length(self.xFloat));
                end
                if ~isempty(self.xDrifter)
                    self.netcdfTool.InitializeDrifterStorage(length(self.xDrifter));
                end
                if ~isempty(self.tracers)
                    for iTracer = 1:length(self.tracers)
                        self.netcdfTool.InitializeTracerStorageWithName(self.tracerNames{iTracer});
                    end
                end
            else
                if isempty(self.netcdfTool.ncid)
                    self.netcdfTool.open();
                end
            end

            self.incrementsWrittenToFile = 0;
            stepsPerOutput = round(self.outputInterval/self.integrator.stepSize);
            firstIncrement = round((self.t0-self.t)/self.integrator.stepSize);
            self.outputIncrements = firstIncrement:stepsPerOutput:self.nIncrements;
%             outputTimes = self.t + self.outputIncrements*self.integrator.stepSize;

            % Save the initial conditions
            self.WriteTimeStepToNetCDFFile(0,self.t);         
        end


        function WriteTimeStepToNetCDFFile(self, integratorIncrement)
            if ~isempty(self.outputFile) && self.timeIndex <= length(self.outputIncrements)
                if self.outputIncrements(self.timeIndex) == integratorIncrement
                    self.netcdfTool.WriteTimeAtIndex(self.timeIndex,self.t);

                    if self.linearDynamics == 0
                        self.netcdfTool.WriteAmplitudeCoefficientsAtIndex(self.timeIndex);
                        self.netcdfTool.WriteEnergeticsAtIndex(self.timeIndex);
                        self.netcdfTool.WriteEnergeticsKJAtIndex(self.timeIndex);
                    end

                    if ~isempty(self.xFloat)
                        self.netcdfTool.WriteFloatPositionsAtIndex(self.timeIndex,self.xFloat,self.yFloat,self.zFloat);
                    end
                    if ~isempty(self.xDrifter)
                        self.netcdfTool.WriteDrifterPositionsAtIndex(self.timeIndex,self.xDrifter,self.yDrifter,self.zDrifter);
                    end
                    if ~isempty(self.tracers)
                        for iTracer = 1:length(self.tracers)
                            self.netcdfTool.WriteTracerWithNameAtIndex(self.timeIndex,self.tracers{iTracer},self.tracerNames{iTracer});
                        end
                    end

                    self.netcdfTool.sync();
                    self.timeIndex = self.timeIndex + 1;
                    self.incrementsWrittenToFile = self.incrementsWrittenToFile + 1;
                end
            end
        end

        function CloseNetCDFFile(self)
            if ~isempty(self.outputFile)
                fprintf('Ending simulation. Wrote %d time points to file\n',self.incrementsWrittenToFile);
                self.netcdfTool.close();
            end
        end

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Integration
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
            if self.linearDynamics == 0
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
            if self.linearDynamics == 0
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