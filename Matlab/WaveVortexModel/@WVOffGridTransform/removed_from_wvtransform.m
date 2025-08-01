        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %
        % Add and remove off-grid internal waves
        %
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        fillOutWaveSpectrum(self,maxTimeGap)
        
        function removeAllExternalWaves(self)
            % remove all external (non-gridded) waves
            %
            % - Topic: External (non-gridded) modes
            self.offgridModes.removeAllExternalWaves();
        end
        
        function omega = setExternalWavesWithWavenumbers(self, k, l, j, phi, A, norm)
            % set external (non-gridded) waves with a given wavenumber
            %
            % - Topic: External (non-gridded) modes
            omega = self.offgridModes.setExternalWavesWithWavenumbers(k, l, j, phi, A, norm);
        end

        function omega = addExternalWavesWithWavenumbers(self, k, l, j, phi, A, norm)
            % add external (non-gridded) waves with a given wavenumber
            %
            % - Topic: External (non-gridded) modes
            omega = self.offgridModes.addExternalWavesWithWavenumbers(k, l, j, phi, A, norm);
        end
        
        function k = setExternalWavesWithFrequencies(self, omega, alpha, j, phi, A, norm)
            % set external (non-gridded) waves with a given frequency
            %
            % - Topic: External (non-gridded) modes
            k = self.offgridModes.setExternalWavesWithFrequencies(omega, alpha, j, phi, A, norm);
        end
        
        function k = addExternalWavesWithFrequencies(self, omega, alpha, j, phi, A, norm)
            % set external (non-gridded) waves with a given wavenumber
            %
            % - Topic: External (non-gridded) modes
            k = self.offgridModes.addExternalWavesWithFrequencies(omega, alpha, j, phi, A, norm);
        end
        
        function [varargout] = externalVariableFieldsAtTime(self,t,varargin)
            % Returns the external wave modes at the grid points.
            %
            % - Topic: External (non-gridded) modes
            varargout = cell(size(varargin));
            [varargout{:}] = self.offgridModes.externalVariablesAtTimePosition(t,reshape(self.X,[],1),reshape(self.Y,[],1), reshape(self.Z,[],1), varargin{:});
            for iArg=1:length(varargout)
                varargout{iArg} = reshape(varargout{iArg},self.Nx,self.Ny,self.Nz);
            end
        end

        function [varargout] = externalVariablesAtTimePosition(self,t,x,y,z,varargin)
        % Returns the external wave modes at the grid points.
        %
        % - Topic: External (non-gridded) modes
        varargout = cell(size(varargin));
        [varargout{:}] = self.offgridModes.externalVariablesAtTimePosition(t,x,y,z, varargin{:});
        end