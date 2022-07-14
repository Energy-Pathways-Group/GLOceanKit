function [zIsopycnal, rhoIsopycnal] = PlaceParticlesOnIsopycnal(self,x,y,z,varargin)
            % MAS 1/10/18 - added intext ('int' or 'both') to give option of using int vs. int+ext fields for rho_prime
            % Also added rhoIsopycnal as output.
            % Any floats with the same value of z will be moved to the same
            % isopycnal.
            %
            % interpolation should probably be 'spline'.
            % tolerance is in meters, 1e-8 should be fine.

            if mod(length(varargin),2) ~= 0
                error('Arguments must be given as name/value pairs.');
            end
            interpolationMethod = 'spline';
            tolerance = 1e-8;
            maxIterations = 200; % max number of iterations to attempt to converge
            useModes = 'all';
            shouldShowDiagnostics = 0;
            for iArg = 1:2:length(varargin)
                if strcmp(varargin{iArg}, 'interpolationMethod')
                    interpolationMethod = varargin{iArg+1};
                elseif strcmp(varargin{iArg}, 'tolerance')
                    tolerance = varargin{iArg+1};
                elseif strcmp(varargin{iArg}, 'maxIterations')
                    maxIterations = varargin{iArg+1};
                elseif strcmp(varargin{iArg}, 'shouldShowDiagnostics')
                    shouldShowDiagnostics = varargin{iArg+1};
                elseif strcmp(varargin{iArg}, 'useModes')
                    if strcmp(varargin{iArg+1}, 'all') || strcmp(varargin{iArg+1}, 'internalOnly') || strcmp(varargin{iArg+1}, 'externalOnly')
                        useModes = varargin{iArg+1};
                    else
                        error('Invalid option for initializeModes');
                    end
                else
                    error('Invalid argument');
                end
            end
            
            t = 0;            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Wrap the particle positions, as necessary
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % (x,y) are periodic for the gridded solution
            x_tilde = mod(x,self.Lx);
            y_tilde = mod(y,self.Ly);
            
            dx = self.x(2)-self.x(1);
            x_index = floor(x_tilde/dx);
            dy = self.y(2)-self.y(1);
            y_index = floor(y_tilde/dy);
            
            % Identify the particles along the interpolation boundary
            if strcmp(interpolationMethod,'spline')
                S = 3+1; % cubic spline, plus buffer
            elseif strcmp(interpolationMethod,'linear')
                S = 1+1;
            end
            bp = x_index < S-1 | x_index > self.Nx-S | y_index < S-1 | y_index > self.Ny-S;
            
            % then do the same for particles along the boundary.
            x_tildeS = mod(x+S*dx,self.Lx);
            y_tildeS = mod(y+S*dy,self.Ly);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %
            % Create the gridded internal density field rho and the interpolants
            %
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Rho = self.InternalVariableFieldsAtTime(t,'rho_prime');
            RhoInterp = griddedInterpolant(self.X,self.Y,self.Z, Rho,interpolationMethod);
            if any(bp)
               RhoInterpS = griddedInterpolant(self.X,self.Y,self.Z, circshift(Rho,[S S 0]),interpolationMethod); 
            end
                    
            % Now let's place the floats along an isopycnal.
            zIsopycnal = z;
            zUnique = unique(z);
            
            if shouldShowDiagnostics == 1
                figure('Position',[100 100 700 600]);
            end
            
            % MAS 1/10/18: create scale factor (<=1.0) to avoid overshoot when trying to converge to isopycnal
            fac = 0.3;
            
            rho = zeros(size(zIsopycnal));                      % initialize rho array
            for zLevel = 1:length(zUnique)
                iterations = 0;
                zLevelIndices = (z==zUnique(zLevel));           % 1 if this is a zLevel, 0 if not
                
                nonboundaryIndices = zLevelIndices & ~bp;       % 1 if this is not a boundary index, 0 if it is
                % get rho for all non-boundary float locations
                rho(nonboundaryIndices) = RhoInterp(x_tilde(nonboundaryIndices),y_tilde(nonboundaryIndices),zIsopycnal(nonboundaryIndices));
                % now get rho for all boundary float locations
                if any(bp)
                    boundaryIndices = zLevelIndices & bp;
                    rho(boundaryIndices) = RhoInterpS(x_tildeS(boundaryIndices),y_tildeS(boundaryIndices),zIsopycnal(boundaryIndices));
                end
                
                % MAS 1/10/18: now make rho = gridded rho_prime + rho_bar (+ external rho_prime)
                if strcmp(useModes,'internalOnly')==1
                  rho(zLevelIndices) = rho(zLevelIndices) + self.RhoBarAtDepth(zIsopycnal(zLevelIndices));
                  elseif strcmp(useModes,'externalOnly')==1
                  rho(zLevelIndices) = self.RhoBarAtDepth(zIsopycnal(zLevelIndices)) + self.externalVariablesAtTimePosition(t,x(zLevelIndices),y(zLevelIndices),zIsopycnal(zLevelIndices),'rho_prime');
                elseif strcmp(useModes,'all')==1
                  rho(zLevelIndices) = rho(zLevelIndices) + self.RhoBarAtDepth(zIsopycnal(zLevelIndices)) + self.externalVariablesAtTimePosition(t,x(zLevelIndices),y(zLevelIndices),zIsopycnal(zLevelIndices),'rho_prime');
                end
                
                dRho = rho(zLevelIndices) - mean(rho(zLevelIndices));
                dz = dRho * 9.81./(self.N2AtDepth(zIsopycnal(zLevelIndices))*self.rho0);
                zIsopycnal(zLevelIndices) = zIsopycnal(zLevelIndices)+dz;
                
                while( max(abs(dz)) > tolerance && iterations < maxIterations )
                    rho(nonboundaryIndices) = RhoInterp(x_tilde(nonboundaryIndices),y_tilde(nonboundaryIndices),zIsopycnal(nonboundaryIndices));
                   if any(bp)
                        rho(boundaryIndices) = RhoInterpS(x_tildeS(boundaryIndices),y_tildeS(boundaryIndices),zIsopycnal(boundaryIndices));
                    end

                    % MAS 1/10/18: now make rho = gridded rho_prime + rho_bar (+ external rho_prime)
                    if strcmp(useModes,'internalOnly')==1
                        rho(zLevelIndices) = rho(zLevelIndices) + self.RhoBarAtDepth(zIsopycnal(zLevelIndices));
                    elseif strcmp(useModes,'externalOnly')==1
                        rho(zLevelIndices) = self.RhoBarAtDepth(zIsopycnal(zLevelIndices)) + self.externalVariablesAtTimePosition(t,x(zLevelIndices),y(zLevelIndices),zIsopycnal(zLevelIndices),'rho_prime');
                    elseif strcmp(useModes,'all')==1
                        rho(zLevelIndices) = rho(zLevelIndices) + self.RhoBarAtDepth(zIsopycnal(zLevelIndices)) + self.externalVariablesAtTimePosition(t,x(zLevelIndices),y(zLevelIndices),zIsopycnal(zLevelIndices),'rho_prime');
                    end

                    dRho = rho(zLevelIndices) - mean(rho(zLevelIndices));
                    dz = dRho * 9.81./(self.N2AtDepth(zIsopycnal(zLevelIndices))*self.rho0);
                    zIsopycnal(zLevelIndices) = zIsopycnal(zLevelIndices)+fac*dz;
                    iterations = iterations + 1;
                    if shouldShowDiagnostics == 1
                        fprintf('  PlaceParticlesOnIsopycnal: Iteration = %d. Mean dz=%3d m of isopycnal at z=%.1f m\n',iterations,mean(abs(dz)),mean(z(zLevelIndices)));
                        if iterations==1
                            subplot(2,1,1)
                            plot(zIsopycnal(zLevelIndices),'ro')
                            hold on
                            subplot(2,1,2)
                            plot(rho(zLevelIndices),'ro')
                            hold on
                        else
                            subplot(2,1,1)
                            plot(zIsopycnal(zLevelIndices),'bo')
                            subplot(2,1,2)
                            plot(rho(zLevelIndices),'bo')
                        end
                    end
                end
                % Do this one more time now that out of 'while' loop to get final rho - we'll pass this back along with zIsopycnal
                rho(nonboundaryIndices) = RhoInterp(x_tilde(nonboundaryIndices),y_tilde(nonboundaryIndices),zIsopycnal(nonboundaryIndices));
                if any(bp)
                  rho(boundaryIndices) = RhoInterpS(x_tildeS(boundaryIndices),y_tildeS(boundaryIndices),zIsopycnal(boundaryIndices));
                end

                % MAS 1/10/18: now make rho = gridded rho_prime + rho_bar (+ external rho_prime)
                if strcmp(useModes,'internalOnly')==1
                    rho(zLevelIndices) = rho(zLevelIndices) + self.RhoBarAtDepth(zIsopycnal(zLevelIndices));
                elseif strcmp(useModes,'externalOnly')==1
                    rho(zLevelIndices) = self.RhoBarAtDepth(zIsopycnal(zLevelIndices)) + self.externalVariablesAtTimePosition(t,x(zLevelIndices),y(zLevelIndices),zIsopycnal(zLevelIndices),'rho_prime');
                elseif strcmp(useModes,'all')==1
                  rho(zLevelIndices) = rho(zLevelIndices) + self.RhoBarAtDepth(zIsopycnal(zLevelIndices)) + self.externalVariablesAtTimePosition(t,x(zLevelIndices),y(zLevelIndices),zIsopycnal(zLevelIndices),'rho_prime');
                end

                if shouldShowDiagnostics == 1
                    subplot(2,1,1)
                    plot(zIsopycnal(zLevelIndices),'m*')
                    ylabel('float depth (m)')
                    subplot(2,1,2)
                    plot(rho(zLevelIndices),'m*')
                    xlabel('Float #')
                    ylabel('float \rho (kg/m^3)')
                end
                
                fprintf('Num iterations=%d. All floats within %.2g m of isopycnal at z=%.1f meters\n',iterations,max(abs(dz)),mean(z(zLevelIndices)) )
            end
            rhoIsopycnal = rho;
        end
