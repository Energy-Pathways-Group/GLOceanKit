classdef WVArrayIntegrator < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        stepSize

        initialTime
        currentTime
        finalTime
        tspan

        totalIterations
        
        OutputFcn = []

        fFromTY
        currentY
        
        previousT
        previousY
        F1, F2, F3, F4
    end
    
    methods
        function self = WVArrayIntegrator( f, tspan, y0, dt, options ) 
            arguments
                f function_handle
                tspan double
                y0 cell
                dt double
                options.OutputFcn
            end
            if length(tspan) < 2
                error('tspan must contain, at minimum, a initial and final time.');
            end
            self.tspan = tspan;
            self.initialTime = tspan(1);
            self.finalTime = tspan(end);

            if isfield(options,"OutputFcn")
                self.OutputFcn = options.OutputFcn;
            end

            f0 = f(self.initialTime,y0);
            for i=1:length(y0)
                if ~isequal(size(y0{i}),size(f0{i}))
                    error('Inconsistent sizes of Y0 and f(t0,y0).');
                end
            end
            
            self.stepSize = dt;
            self.totalIterations = 0;
            self.fFromTY = f;
            self.currentY = y0;  
            self.currentTime = self.initialTime;

            self.previousT = self.currentTime;
            self.previousY = self.currentY;

            self.integrateToTime(self.finalTime);
        end
            
        function integrateToTime(self, time )
            if ~isempty(self.OutputFcn)
                self.OutputFcn(self.currentTime,self.currentY,'init');
            end
            while self.currentTime < time                
                self.currentY = self.stepForward(self.currentY,self.currentTime,self.stepSize);
                self.totalIterations = self.totalIterations + 1;
                self.currentTime = self.initialTime + self.totalIterations * self.stepSize;
                if ~isempty(self.OutputFcn)
                    t_output = self.tspan(self.tspan > self.previousT & self.tspan <= self.currentTime);
                    for iOutput = 1:length(t_output)
                        self.OutputFcn(t_output(iOutput),self.valueAtTime(t_output(iOutput)),'');
                    end
                end
            end
            if ~isempty(self.OutputFcn)
                self.OutputFcn(self.finalTime,self.valueAtTime(self.finalTime),'done');
            end
        end
   
        function yo = stepForward(self,yi,t,dt)
            self.previousY = yi;
            self.previousT = t;
            self.F1 = feval(self.fFromTY,t,yi);
            
            y2 = cell(size(yi));
            for i=1:length(y2)
                y2{i} = yi{i}+0.5*dt*self.F1{i};
            end
            self.F2 = feval(self.fFromTY,t+0.5*dt,y2);
            
            y3 = cell(size(yi));
            for i=1:length(y3)
                y3{i} = yi{i}+0.5*dt*self.F2{i};
            end
            self.F3 = feval(self.fFromTY,t+0.5*dt,y3);
            
            y4 = cell(size(yi));
            for i=1:length(y4)
                y4{i} = yi{i}+dt*self.F3{i};
            end
            self.F4 = feval(self.fFromTY,t+dt,y4);
            
            yo = cell(size(yi));
            for i=1:length(yo)
                yo{i} = yi{i} + (dt/6)*(self.F1{i} + 2*self.F2{i} + 2*self.F3{i} + self.F4{i});
            end
        end

        function yo = valueAtTime(self,t)
            % Hermite interpolation
            theta = (t - self.previousT)/self.stepSize;
            if theta < -eps || theta > 1+eps
                error("invalid time for interpolation");
            end
            alpha_2 = 3*theta*theta - 2*theta*theta*theta;
            alpha_1 = 1 - alpha_2;
            alpha_3 = self.stepSize*(theta - 2*theta*theta + theta*theta*theta);
            alpha_4 = self.stepSize*(theta*theta*theta - theta*theta);

            yo = cell(size(self.currentY));
            for i = 1:numel(yo)
                yo{i} = alpha_1*self.previousY{i} + alpha_2*self.currentY{i} + alpha_3*self.F1{i} + alpha_4*self.F4{i};
            end
        end
        
    end
    
end

