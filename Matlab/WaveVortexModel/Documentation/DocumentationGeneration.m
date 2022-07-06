mc = meta.class.fromName('WaveVortexModel');

for i=1:length(mc.MethodList)
    if strcmp(mc.MethodList(i).DefiningClass.Name,'handle')
        continue;
    end

    if ~strcmp(mc.MethodList(i).Access,'public') || (mc.MethodList(i).Hidden == true)
        continue;
    end
        fprintf('%d: %s\n',i,mc.MethodList(i).Name)
end

for i=1:length(mc.PropertyList)
    if strcmp(mc.PropertyList(i).DefiningClass.Name,'handle')
        continue;
    end

    if ~strcmp(mc.PropertyList(i).GetAccess,'public')
        continue;
    end
        fprintf('%d: %s\n',i,mc.PropertyList(i).Name)
end