% Load our giant list of possible paths.
K = 4;
dof = 4;
totalPermutations = 1000;
shouldSaveFigures = 0;


figure('Units', 'points', 'Position', [50 50 figure_width_1col 175*scaleFactor])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

for SiteNumber = 1:2

    load(sprintf('../observations/smoothedGriddedRho%dDrifters.mat',SiteNumber));
    load(sprintf('Rho%dDrifterSplineFits%d_K%d_dof%d.mat',SiteNumber,totalPermutations,K,dof));
    
    scaleFactor = 1;
    LoadFigureDefaults
    
    f0 = 2 * 7.2921E-5 * sin( lat0*pi/180. );
    x = x(:,1:(end-1));
    y = y(:,1:(end-1));
    nDrifters = size(x,2);
    nT = size(x,1);
    D = FiniteDifferenceMatrix(1,t,1,1,2);
    mx = mean(x,2);
    my = mean(y,2);
    dmxdt = D*mx;
    dmydt = D*my;
    dxdt = D*x;
    dydt = D*y;
    
    if SiteNumber == 1
        iModel = 3;
    elseif SiteNumber == 2
        iModel = 5;
    end
    
    [~,mostLikelyIndices] = sort(bootstraps{iModel}.jointlikelihood,'descend');
    indexBest = mostLikelyIndices(1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %
    % Decomposition figures
    %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    p = bootstraps{iModel};
    u0 = p.u0(:,indexBest);
    v0 = p.v0(:,indexBest);
    sigma_n = p.sigma_n(:,indexBest);
    sigma_s = p.sigma_s(:,indexBest);
    zeta = p.zeta(:,indexBest);
    delta = p.delta(:,indexBest);
    
    u_meso = u0 + 0.5*(sigma_n+delta).*x + 0.5*(sigma_s-zeta).*y;
    v_meso = v0 + 0.5*(sigma_s+zeta).*x + 0.5*(delta-sigma_n).*y;
    x_meso = x(1,:) + cumtrapz(t,u_meso);
    y_meso = y(1,:) + cumtrapz(t,v_meso);
    
    u_bg = dmxdt - mean(u_meso,2);
    v_bg = dmydt - mean(v_meso,2);
    x_bg = cumtrapz(t,u_bg);
    y_bg = cumtrapz(t,v_bg);
    
    u_sm = dxdt - u_meso - u_bg;
    v_sm = dydt - v_meso - v_bg;
    x_sm = cumtrapz(t,u_sm);
    y_sm = cumtrapz(t,v_sm);
    
    
    q = (x-mx);
    r = (y-my);
    dqdt = D*q;
    drdt = D*r;
    dqdt_meso = 0.5*(sigma_n+delta).*q + 0.5*(sigma_s-zeta).*r;
    drdt_meso = 0.5*(sigma_s+zeta).*q + 0.5*(delta-sigma_n).*r;
    
    q_meso = q(1,:) + cumtrapz(t,dqdt_meso);
    r_meso = r(1,:) + cumtrapz(t,drdt_meso);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Quick aside
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    measuredTotalVariance = sum(dqdt.^2 + drdt.^2,2);
    modelMesoscaleVariance = sum(dqdt_meso.^2 + drdt_meso.^2,2);
    modelSubmesoscaleVariance = sum(u_sm.^2 + v_sm.^2,2);
    modelTotalVariance = sum(dqdt_meso.^2 + drdt_meso.^2+u_sm.^2 + v_sm.^2,2);
    EnergyCrossTerms = (measuredTotalVariance - modelTotalVariance)./measuredTotalVariance;
    % figure
    % plot(t/86400,EnergyCrossTerms) %, hold on, plot(t/86400,actualTotalVariance);
    
    
    cv_ms = dqdt_meso + sqrt(-1)*drdt_meso;
    cv_sm = u_sm + sqrt(-1)*v_sm;
    
    [psi,lambda]=sleptap(size(cv_ms,1));
    dt = t(2)-t(1);
    [f,sxx,syy,sxy]=mspec(dt,cv_ms,cv_sm,psi,'cyclic');
    gamma=frac(abs(sxy).^2,sxx.*syy);
    
    
    plot(f*86400,mean(gamma,2),'LineWidth',2), hold on
    xlim([0 12])
    ylim([0 1])
    xlabel('frequency (cycles/day)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
    ylabel('coherence', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
    p1 = gca;
    p1.FontSize = figure_axis_tick_size;
end

legend('Site 1', 'Site 2')

if shouldSaveFigures == 1
    print('CoherenceSite1Site2.eps','-depsc2');
end
