function [flag,xi,yi] = doesIntersect(seg1,seg2)
% seg1 = [x1 y1; x2 y2];
% seg3 = [x3 y3; x4 y4];
% https://en.wikipedia.org/wiki/Line–line_intersection
x1 = seg1(1,1); y1 = seg1(1,2);
x2 = seg1(2,1); y2 = seg1(2,2);
x3 = seg2(1,1); y3 = seg2(1,2);
x4 = seg2(2,1); y4 = seg2(2,2);

flag = 0;

denom = (x1-x2)*(y3-y4)-(y1-y2)*(x3-x4);
if denom == 0    
    return;
end

a = (x1-x3)*(y3-y4) - (y1-y3)*(x3-x4);
b = (x1-x2)*(y1-y3) - (y1-y2)*(x1-x3);

t = a/denom;
u = -b/denom;

flag = t >= 0 && t <= 1 && u >= 0 && u <= 1;

if flag == 1
   xi = x1 + t*(x2-x1);
   yi = y1 + t*(y2-y1);
else
    xi = NaN;
    yi = NaN;
end

end