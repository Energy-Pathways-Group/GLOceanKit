function dLevels = DensityLevelForCDF(X,Y,density, pctTarget)
% Now let's identify the cumulative values enclosed by the different
% contours.
nLevels = 25;
level = zeros(nLevels,1);
pctEnclosed = zeros(nLevels,1);
sigma_n_axis = X(1,:).';
sigma_s_axis = Y(:,1);
M = contourc(sigma_n_axis,sigma_s_axis,density,nLevels);
i = 1;
iLevel = 1;
while (i < size(M,2))
    level(iLevel) = M(1,i);
    n = M(2,i);
    
    % the last point is redudant, and polyshape doesn't like that.
%     cont = polyshape(M(1,(i+1):(i+n-1)),M(2,(i+1):(i+n-1)));
%     mask = isinterior(cont,reshape(X,[],1),reshape(Y,[],1));
    
    cont.Vertices = cat(2,M(1,(i+1):(i+n-1)).',M(2,(i+1):(i+n-1)).');
    cont.dx = cont.Vertices(2:end,1) - cont.Vertices(1:(end-1),1);
    cont.dy = cont.Vertices(2:end,2) - cont.Vertices(1:(end-1),2);
    mask = isInterior(cont,reshape(X,[],1),reshape(Y,[],1));
    
    mask = reshape(mask,size(X));
    
    pctEnclosed(iLevel) = trapz(sigma_s_axis,trapz(sigma_n_axis,density.*mask));
    
    i = i+n+1;
    iLevel = iLevel + 1;
end

dLevels = interp1(pctEnclosed,level,pctTarget);
end

function isLeft = isLeft(polygon,i,x,y)
% using the notation above,
% -x3x4*y2y3 + x2x3*y3y4
isLeft = (polygon.dx(i)*(y-polygon.Vertices(i,2)) - (x-polygon.Vertices(i,1))*polygon.dy(i) );
end

function isinterior = isInterior(polygon, x, y)
% In my tests this was a 6-7x speedup over matlabs
% implementation.
windingNumber = zeros(size(x));
for i=1:(length(polygon.Vertices)-1)
    isless = polygon.Vertices(i,2) <= y;
    upwardCrossing = isless & polygon.Vertices(i+1,2) > y;
    if any(upwardCrossing)
        windingNumber(upwardCrossing) = windingNumber(upwardCrossing) + (isLeft(polygon,i,x(upwardCrossing),y(upwardCrossing)) > 0);
    end
    downwardCrossing = ~isless & polygon.Vertices(i+1,2) <= y;
    if any(downwardCrossing)
        windingNumber(downwardCrossing) = windingNumber(downwardCrossing) - (isLeft(polygon,i,x(downwardCrossing),y(downwardCrossing)) < 0);
    end
end
isinterior = abs(windingNumber) > 0;
end