% create a box

L = 5;
xbounds = [-Inf Inf];
ybounds = [-Inf Inf]; %[-1 1];

xbounds = 2*[-L L];
ybounds = 2*[-L L]; %[-1 1];

% outer box
xv = [min(xbounds) min(xbounds) max(xbounds) max(xbounds)];
yv = [max(ybounds) min(ybounds) min(ybounds) max(ybounds)];

pgon = polyshape(xv,yv);

xbounds = [-L L];
ybounds = [-L L];

% inner box
xv_i = [min(xbounds) max(xbounds) max(xbounds) min(xbounds)];
yv_i = [max(ybounds) max(ybounds) min(ybounds) min(ybounds)];

pgon = addboundary(pgon,xv_i,yv_i);

% in = inpolygon(-100,2,xv,yv)

xv = cat(2,xv,NaN, xv_i);
yv = cat(2,yv,NaN, yv_i);


theta = linspace(0,2*pi-2*pi/30,30);
pgon = addboundary(pgon,cos(theta),sin(theta));

% xv = cat(2,xv,NaN, cos(theta));
% yv = cat(2,yv,NaN, sin(theta));

% pgon = polyshape(xv,yv);
figure, plot(pgon)

isinterior(pgon,0,3)
isinterior(pgon,0,0)