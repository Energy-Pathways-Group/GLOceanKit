mc = meta.class.fromName('WaveVortexModel');

for i=1:length(mc.MethodList)
    if strcmp(mc.MethodList(i).DefiningClass.Name,'handle')
        continue;
    end

    if ~strcmp(mc.MethodList(i).Access,'public')
        continue;
    end
        fprintf('%d: %s\n',i,mc.MethodList(i).Name)
end