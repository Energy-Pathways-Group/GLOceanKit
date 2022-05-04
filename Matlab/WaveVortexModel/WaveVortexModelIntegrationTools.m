classdef WaveVortexModelIntegrationTools < handle

    properties
        wvm             % WaveVortexModel
        t=0             % current model time (in seconds)
    
        outputFile
        netcdfTool      % WaveVortexModelNetCDFTools instance---empty indicates no file output
        outputInterval  % model output interval (seconds)
        t0=0            % initial output time
        iTime=1         % current index of outputTimes written to file.

        integrator      % Array integrator
    end

    methods

        function self = IntegrationToolForModel(self,waveVortexModel)
            self.wvm = waveVortexModel;
        end

        function self = IntegrationToolForModelRestart(self,existingModelOutput,restartIndex,shouldDoubleResolution)
        %IntegratorForModelRestart Initializes a new instance of the
        %   WaveVortexModel using existing output.
        %   existingModelOutput     path to existing netcdf file
        %   restartIndex            time point to start from. Inf will
        %                           restart from the last time point (default).
        %   shouldDoubleResolution  double the model resolution.


        end

        function self = SetNetCDFFileForModelOutput(self,outputFile,outputInterval)
            % outputInterval    Desired output interval
            self.outputFile = outputFile;
            self.outputInterval = outputInterval;
        end


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
                        fprintf('\ttime step %d of %d. Estimated finish time %s (%s from now)\n', i, nIncrements, datestr(datetime('now')+timeRemaining), datestr(timeRemaining, 'HH:MM:SS')) ;
                        self.wvm.summarizeEnergyContent();
                    end
                end

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
                        self.iTime = self.iTime + 1;
                        incrementsWrittenToFile = incrementsWrittenToFile + 1;
                    end
                end
            end

            if ~isempty(self.outputFile)
                fprintf('Ending simulation. Wrote %d time points to file\n',incrementsWrittenToFile);
            end
% closing and re-opening is leading to an error:
% The NetCDF library encountered an error during execution of 'putVaraDouble' function - 'HDF error (NC_EHDFERR)'
% Lots of people started having this issue, and there is no explanation for
% what is going wrong or what the work around is.
%             if ~isempty(self.outputFile)
%                 self.netcdfTool.close();
%             end
        end

        function self = SetupNetCDFToolsForOutputFile(self)
            if isempty(self.outputFile)
                % user didn't request output, so move on
                return;
            end

            if isempty(self.netcdfTool)
                self.netcdfTool = WaveVortexModelNetCDFTools(self.outputFile);
                self.netcdfTool.CreateNetCDFFileFromModel(self.wvm,Inf,'double');
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
            self.integrator = ArrayIntegrator(@(t,y0) self.wvm.NonlinearFluxAtTimeArray(t,y0),{self.wvm.Ap,self.wvm.Am,self.wvm.A0},deltaT);
            self.integrator.currentTime = initialTime;
        end
    end
end