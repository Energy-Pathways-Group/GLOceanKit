classdef WVArrayIntegrator < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        stepSize
        currentTime
        totalIterations
        
        fFromTY
        currentY
        
        F
    end
    
    methods
        function self = WVArrayIntegrator( f, y0, dt, options ) 
            arguments
                f function_handle
                y0 cell
                dt double
                options.currentTime = 0
            end
f0 = feval(f,options.currentTime,y0);
%             try
%                 f0 = feval(f,0,y0);
%             catch theError
%                 msg = ['Unable to evaluate the ODEFUN at t0,y0. ',theError];
%                 error(msg);
%             end
            
            for i=1:length(y0)
                if ~isequal(size(y0{i}),size(f0{i}))
                    error('Inconsistent sizes of Y0 and f(t0,y0).');
                end
            end
            
            self.currentTime = 0;
            self.stepSize = dt;
            self.totalIterations = 0;
            self.fFromTY = f;
            self.currentY = y0;  
            self.currentTime = options.currentTime;
        end
        
        function [y, t] = IntegrateAlongDimension(self,time)
            y = zeros(self.nReps,self.nDims,length(time));
            t = zeros(size(time));
            for i=1:length(time)
                p = self.StepForwardToTime(time(i));
                y(:,:,i) = p;
                t(i) = self.currentTime;
            end
        end
        
        function y = StepForwardToTime(self, time )
            while self.currentTime < time                
                self.currentY = self.StepForward(self.currentY,self.currentTime,self.stepSize);
                self.currentTime = self.currentTime + self.stepSize;
                self.totalIterations = self.totalIterations + 1;
            end
            
            y = self.currentY;
        end
        
        function y = IncrementForward(self)
            self.currentY = self.StepForward(self.currentY,self.currentTime,self.stepSize);
            self.currentTime = self.currentTime + self.stepSize;
            self.totalIterations = self.totalIterations + 1;
            
            y = self.currentY;
        end
        
        function yo = StepForward(self,yi,t,dt)
            F1 = feval(self.fFromTY,t,yi);
            
            y2 = cell(size(yi));
            for i=1:length(y2)
                y2{i} = yi{i}+0.5*dt*F1{i};
            end
            F2 = feval(self.fFromTY,t+0.5*dt,y2);
            
            y3 = cell(size(yi));
            for i=1:length(y3)
                y3{i} = yi{i}+0.5*dt*F2{i};
            end
            F3 = feval(self.fFromTY,t+0.5*dt,y3);
            
            y4 = cell(size(yi));
            for i=1:length(y4)
                y4{i} = yi{i}+dt*F3{i};
            end
            F4 = feval(self.fFromTY,t+dt,y4);
            
            yo = cell(size(yi));
            for i=1:length(yo)
                yo{i} = yi{i} + (dt/6)*(F1{i} + 2*F2{i} + 2*F3{i} + F4{i});
            end
        end
        
    end
    
end

