classdef WaveVortexModelNetCDFFile < NetCDFFile

    properties
        wvm     % waveVortexModel instance
        t       % model time of the wvm coefficients
    end

    methods
        function self = WaveVortexModelNetCDFFile(path,overwriteExisting)
            self@NetCDFFile(path,overwriteExisting);

        end

        function CreateNetCDFFileFromModel(self,waveVortexModel,Nt,precision)
            self.wvm = waveVortexModel;

            self.addDimension(self.wvm.x,'x',containers.Map({'units'},{'m'}));
            self.addDimension(self.wvm.y,'y',containers.Map({'units'},{'m'}));
            self.addDimension(self.wvm.z,'z',containers.Map({'units'},{'m'}));
            self.addDimension(self.wvm.k,'k',containers.Map({'units'},{'radians/m'}));
            self.addDimension(self.wvm.l,'l',containers.Map({'units'},{'radians/m'}));
            self.addDimension(self.wvm.j,'j',containers.Map({'units'},{'mode number'}));
            self.addDimension([],'t',containers.Map({'units'},{'s'}),Nt);

            self.addVariable(self.wvm.IMA0,'IMA0',{'k','l','j'});
            self.addVariable(self.wvm.IMAp,'IMAp',{'k','l','j'});
            self.addVariable(self.wvm.IMAm,'IMAm',{'k','l','j'});
            self.addVariable(self.wvm.EMA0,'EMA0',{'k','l','j'});
            self.addVariable(self.wvm.EMAp,'EMAp',{'k','l','j'});
            self.addVariable(self.wvm.EMAm,'EMAm',{'k','l','j'});

            %%% !!! Need to make a custom setter for attributes to actually
            %%% write this to file! And probably do the same for variables,
            %%% etc. Or, actually, addAttribute is probably more
            %%% appropriate
            CreationDate = datestr(datetime('now'));
            self.attributes('latitude') = self.wvm.latitude;
            self.attributes('t0') = self.wvm.t0;
            self.attributes('rho0') = self.wvm.rho0;
            self.attributes('Model') = 'Created from WaveVortexModel.m written by Jeffrey J. Early.';
            self.attributes('ModelVersion') = self.wvm.version;
            self.attributes('CreationDate') = CreationDate;

            if isa(self.wvm,'WaveVortexModelConstantStratification')
                self.attributes('stratification') = 'constant';
                self.attributes('N0') = self.wvm.N0;
            elseif isa(self.wvm,'WaveVortexModelHydrostatic')
                self.attributes('stratification') = 'custom-hydrostatic';

                self.addVariable(self.wvm.rhobar,'rhobar',{'z'});
                self.addVariable(self.wvm.N2,'N2',{'z'});
                self.addVariable(self.wvm.dLnN2,'dLnN2',{'z'});

                self.addVariable(self.wvm.PFinv,'PFinv',{'z','j'});
                self.addVariable(self.wvm.QGinv,'QGinv',{'z','j'});
                self.addVariable(self.wvm.PF,'PF',{'j','z'});
                self.addVariable(self.wvm.QG,'QG',{'j','z'});
                self.addVariable(self.wvm.h,'h',{'j'});
                self.addVariable(self.wvm.P,'P',{'j'});
                self.addVariable(self.wvm.Q,'Q',{'j'});

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