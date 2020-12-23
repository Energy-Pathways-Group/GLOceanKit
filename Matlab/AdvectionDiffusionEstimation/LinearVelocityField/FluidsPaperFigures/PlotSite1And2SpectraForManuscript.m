% AMS figure widths, given in picas, converted to points (1 pica=12 points)
scaleFactor = 1;
LoadFigureDefaults

figure('Units', 'points', 'Position', [50 50 figure_width_1col+22 250*scaleFactor])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

for Site=1:2
% Choose which drifters to analyze.
    if Site == 1
        dof=4;
        iModel=3;
        ylimit = [2e-3 1.5e3];
    else
        dof=4;
        iModel=5;
        ylimit = [2e-2 5e4];
    end
    
    load(sprintf('smoothedGriddedRho%dDrifters.mat',Site));
    % The last drifter at both sites is only partial time series.
    x = x(:,1:(end-1));
    y = y(:,1:(end-1));

    filename = sprintf('BootstrapData/Rho%dDrifterSplineFits1000_dof%d.mat',Site,dof);
    load(filename);
    [~,mostLikelyIndices] = sort(bootstraps{iModel}.jointlikelihood,'descend');
    p = bootstraps{iModel};
    parameterEstimates.u0 = p.u0(:,mostLikelyIndices(1));
    parameterEstimates.v0 = p.v0(:,mostLikelyIndices(1));
    parameterEstimates.u1 = p.u1(:,mostLikelyIndices(1));
    parameterEstimates.v1 = p.v1(:,mostLikelyIndices(1));
    parameterEstimates.sigma_n = p.sigma_n(:,mostLikelyIndices(1));
    parameterEstimates.sigma_s = p.sigma_s(:,mostLikelyIndices(1));
    parameterEstimates.zeta = p.zeta(:,mostLikelyIndices(1));
    parameterEstimates.delta = p.delta(:,mostLikelyIndices(1));
    
    [u_meso,v_meso,u_bg,v_bg,u_sm,v_sm,dmxdt,dmydt] = DecomposeTrajectories(x, y, t, parameterEstimates);
    
    cv_bg = u_bg + sqrt(-1)*v_bg;
    cv_strain = u_meso + sqrt(-1)*v_meso;
    cv_sm = u_sm + sqrt(-1)*v_sm;


    %
    % Compute energy spectrum of the center of mass.
    %
    dt = t(2)-t(1);
    [psi,lambda]=sleptap(size(cv_bg,1),3);
    [~,spp_bg,snn_bg,spn_bg]=mspec(dt, cv_bg,psi);
    [~,spp_strain,snn_strain,spn_strain]=mspec(dt, cv_strain,psi);
    [f,spp_sm,snn_sm,spn_sm]=mspec(dt, cv_sm,psi);

    % convert from radians/s to cycles/day
    f = f*86400/(2*pi);

    dim1 = 0.0; % Setting this above zero allows you to 'dim' the color of the line.

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot tides and inertial
    M2TidalFrequency = 24/12.421;
    S2TidalFrequency = 24/12.00;
    K1TidalFrequency = 24/23.934;
    O1TidalFrequency = 24/25.819;
    
    sp1 = subplot(2,1,Site);
    a = plot(-corfreq(lat0)*24/(2*pi)*ones(1,2),ylimit, 'LineWidth', 1.0*scaleFactor, 'Color', 0.3*[1.0 1.0 1.0] );
    hold on
    b = plot(M2TidalFrequency*ones(1,2),ylimit, 'LineWidth', 1.0*scaleFactor, 'Color', 0.3*[1.0 1.0 1.0]);
    ylog
    xlim([-max(f)/2 max(f)/2])
    ylim(ylimit)
    ylabel('power (m^2/s)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Plot background
    fig_bg = plot([-flip(f,1); f], [flip(snn_bg,1); spp_bg], 'LineWidth', 1.5, 'Color', [dim1 dim1 dim1]);

    dim2 = 0.0;
    fig_strain = plot([-flip(f,1); f], vmean([flip(snn_strain,1); spp_strain],2), 'LineWidth', 1.5, 'Color', [dim2 dim2 1.0]);

    dim3 = 0.0;
    fig_sm = plot([-flip(f,1); f], [flip(vmean(snn_sm,2),1); vmean(spp_sm,2)], 'LineWidth', 1.5, 'Color', [1.0 dim3 1.0]);
    % legend([fig_bg, fig_strain, fig_sm], 'Background', 'Mesoscale', 'Submesoscale', 'Location', 'northeast')

    sp1.YTick = [1e-2 1e0 1e2 1e4];
    sp1.FontSize = figure_axis_tick_size;
    sp1.FontName = figure_font;
    if Site == 1
        sp1.XTickLabel = [];
    else
        xlabel('frequency (cycles per day)', 'FontSize', figure_axis_label_size, 'FontName', figure_font)
    end
end

packfig(2,1)
print('-depsc', 'site1and2_decomposed_spectra.eps')