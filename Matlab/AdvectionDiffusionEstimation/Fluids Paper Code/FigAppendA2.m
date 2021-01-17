load("smoothedGriddedRho1Drifters.mat")
NX = 9;
xl = x(1,1:NX); yl = y(1,1:NX);

delta = 0;  
theta = 30;
zeta = 0;
sigma = 7e-6;
sigma_s = sigma*sind(2*theta);
sigma_n = sigma*cosd(2*theta);

    cmx = mean(xl); cmy = mean(yl);
    cmd = mean(sqrt((xl-cmx).^2+(yl-cmy).^2));
    spc = 9*cmd/20;
    x1 = (-4:4)*spc*cosd(30)+cmx; y1 = (-4:4)*spc*sind(30)+cmy;
    x2 = (-4:4)*spc*cosd(120)+cmx; y2 = (-4:4)*spc*sind(120)+cmy;
    x3=2*(xl-cmx)+cmx; y3=2*(yl-cmy)+cmy;
    x4=.5*(xl-cmx)+cmx; y4=.5*(yl-cmy)+cmy;
    
    a1a=figure; 
    [xq,yq]=meshgrid(-1300:100:1300,-1300:100:1300);
    u_mes = 0.5*(sigma_n+delta).*xq+0.5*(sigma_s-zeta).*yq;
    v_mes = 0.5*(sigma_s+zeta).*xq+0.5*(delta-sigma_n).*yq;
    quiver(xq,yq,u_mes,v_mes);
    
    hold on; scatter(xl-cmx,yl-cmy,'filled','b')
    hold on; scatter(x1-cmx,y1-cmy,'filled','r')
    hold on; scatter(x2-cmx,y2-cmy,'filled','g')
    xlabel('x(m)'); ylabel('y(m)');
    xlim([-1100 1100]); ylim([-1100 1100])
    axis equal
    exportfig(a1a, 'figures/append1a.eps', 'width', 16, 'color', 'cmyk','Fontmode','fixed','FontSize', 14);
    
    a1b=figure;
    [xq,yq]=meshgrid(-2400:200:2400,-2400:200:2400);
    u_mes = 0.5*(sigma_n+delta).*xq+0.5*(sigma_s-zeta).*yq;
    v_mes = 0.5*(sigma_s+zeta).*xq+0.5*(delta-sigma_n).*yq;
    quiver(xq,yq,u_mes,v_mes);
    hold on; scatter(xl-cmx,yl-cmy,'filled','b')
    hold on; scatter(x3-cmx,y3-cmy,'filled','r')
    hold on; scatter(x4-cmx,y4-cmy,'filled','g')
    xlabel('x(m)'); ylabel('y(m)');
    xlim([-2100 2100]); ylim([-2100 2100])
    axis equal
    exportfig(a1b, 'figures/append1b.eps', 'width', 16, 'color', 'cmyk','Fontmode','fixed','FontSize', 14);