classdef WaveVortexModelNetCDFFile < NetCDFFile

    properties
        wvm     % waveVortexModel instance
        t       % model time of the wvm coefficients

        Nx, Ny, Nz
        Nk, Nl, Nj, Nt
        Nkh

        nFloats = 0
        nDrifters = 0
        nTracers = 0
    end

    methods
        function self = WaveVortexModelNetCDFFile(path,varargin)
            if nargin == 0

            elseif isa(varargin{1},'char' ) % overwriteExisting
            self@NetCDFFile(path,varargin{1});
            elseif isa(varargin{2},"double")

            end

        end

        function InitializeFromExistingFile(self)
            InitializeFromExistingFile@NetCDFFile(self);

            requiredVariables = {'x','y','z','j'};
            requiredAttributes = {'latitude','stratification','t0'};
            if ~all(isKey(self.variables,requiredVariables)) || ~all(isKey(self.attributes,requiredAttributes))
                error('This files is missing required variables or attributes to load directly into the WaveVortexModel.')
            end


        end

        function CreateNetCDFFileFromModel(self,waveVortexModel,Nt,precision)
            self.wvm = waveVortexModel;

            self.addDimension('x',self.wvm.x,containers.Map({'units'},{'m'}));
            self.addDimension('y',self.wvm.y,containers.Map({'units'},{'m'}));
            self.addDimension('z',self.wvm.z,containers.Map({'units'},{'m'}));
            self.addDimension('k',self.wvm.k,containers.Map({'units'},{'radians/m'}));
            self.addDimension('l',self.wvm.l,containers.Map({'units'},{'radians/m'}));
            self.addDimension('j',self.wvm.j,containers.Map({'units'},{'mode number'}));
            self.addDimension('t',[],containers.Map({'units'},{'s'}),Nt);

            self.addVariable('IMA0',self.wvm.IMA0,{'k','l','j'});
            self.addVariable('IMAp',self.wvm.IMAp,{'k','l','j'});
            self.addVariable('IMAm',self.wvm.IMAm,{'k','l','j'});
            self.addVariable('EMA0',self.wvm.EMA0,{'k','l','j'});
            self.addVariable('EMAp',self.wvm.EMAp,{'k','l','j'});
            self.addVariable('EMAm',self.wvm.EMAm,{'k','l','j'});

            CreationDate = datestr(datetime('now'));
            self.addAttribute('latitude', self.wvm.latitude);
            self.addAttribute('t0', self.wvm.t0);
            self.addAttribute('rho0',self.wvm.rho0);
            self.addAttribute('Model','Created from WaveVortexModel.m written by Jeffrey J. Early.');
            self.addAttribute('ModelVersion',self.wvm.version);
            self.addAttribute('CreationDate',CreationDate);

            if isa(self.wvm,'WaveVortexModelConstantStratification')
                self.addAttribute('stratification','constant');
                self.addAttribute('N0',self.wvm.N0);
            elseif isa(self.wvm,'WaveVortexModelHydrostatic')
                self.addAttribute('stratification','custom-hydrostatic');

                self.addVariable('rhobar',self.wvm.rhobar,{'z'});
                self.addVariable('N2',self.wvm.N2,{'z'});
                self.addVariable('dLnN2',self.wvm.dLnN2,{'z'});

                self.addVariable('PFinv',self.wvm.PFinv,{'z','j'});
                self.addVariable('QGinv',self.wvm.QGinv,{'z','j'});
                self.addVariable('PF',self.wvm.PF,{'j','z'});
                self.addVariable('QG',self.wvm.QG,{'j','z'});
                self.addVariable('h',self.wvm.h,{'j'});
                self.addVariable('P',self.wvm.P,{'j'});
                self.addVariable('Q',self.wvm.Q,{'j'});

                rhoFunction = self.wvm.rhoFunction;
                N2Function = self.wvm.N2Function;
                dLnN2Function = self.wvm.dLnN2Function;
                save(self.matFilePath,'rhoFunction','N2Function','dLnN2Function','CreationDate');
            else
                error('Not implemented');
            end


        end
    end

end