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
        end
    end

end