theta = linspace(0,2*pi-2*pi/30,30);
poly1 = polyshape(cos(theta),sin(theta));

% poly1 = polyshape([0 0 1 1],[1 0 0 1]);
lineseg = [-1 1; 0 0.5];
lineseg = [-1 1; -0.5 0];
[in,out] = intersect(poly1,lineseg);

figure
plot(poly1)
hold on
% plot(in(:,1),in(:,2),'b',out(:,1),out(:,2),'r')
axis equal

for i=1:length(poly1.Vertices)
    seg1 = lineseg;
    if i == length(poly1.Vertices)
        seg2 = [poly1.Vertices(end,1:2); poly1.Vertices(1,1:2)];
    else
        seg2 = poly1.Vertices(i:(i+1),1:2);
    end
    [flag,xi,yi] = doesIntersect( seg1, seg2 );
    if flag == 1
        plot(seg1(:,1),seg1(:,2),'b','LineWidth', 2)
        plot(seg2(:,1),seg2(:,2),'r','LineWidth', 2)
        scatter(xi,yi,10^2,0*[1 1 1],'filled')
        
        % treat the point of intersection as the origin. This means the
        % interior point is
        px = seg1(2,1) - xi;
        py = seg1(2,2) - yi;
        
        dx = seg2(2,1)-seg2(1,1);
        dy = seg2(2,2)-seg2(1,2);
        
        a = dx*dx-dy*dy;
        b = 2*dx*dy;
        d = dx*dx+dy*dy;
        
        qx = (a*px + b*py)/d + xi;
        qy = (b*px - a*py)/d + yi;
        
        seg3 = [xi yi; qx qy];
        plot(seg3(:,1),seg3(:,2),'g','LineWidth', 2)
    end
end


% reflection
% https://math.stackexchange.com/questions/525082/reflection-across-a-line

