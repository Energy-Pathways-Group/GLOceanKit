% Load our giant list of possible paths.
SiteNumber=2;
totalPermutations = 1000;
load(sprintf('smoothedGriddedRho%dDrifters.mat',SiteNumber));
% The last drifter at both sites is only partial time series.
x = x(:,1:(end-1));
y = y(:,1:(end-1));
shouldFitToSecondMomentOnly = 0;

FVU = @(u_sm,v_sm,u,v) mean(mean(u_sm.^2+v_sm.^2,2))./mean(mean(u.^2+v.^2,2));
FDU = @(u_sm,v_sm,u,v) mean(abs(sum(u_sm+1i*v_sm,1)).^2)/mean(abs(sum(u+1i*v,1)).^2);
kappa = @(u_sm,v_sm) (t(2)-t(1))*mean(abs(sum(u_sm+1i*v_sm,1)).^2)/size(u_sm,1)/4;


for dof = 1:6
    if shouldFitToSecondMomentOnly == 1
        filename = sprintf('BootstrapData/Rho%dDrifterSpline2ndMomFits%d_dof%d.mat',SiteNumber,totalPermutations,dof);
    else
        filename = sprintf('BootstrapData/Rho%dDrifterSplineFits%d_dof%d.mat',SiteNumber,totalPermutations,dof);
    end
    load(filename);
        
    for iModel = 1:length(bootstraps)
        [~,mostLikelyIndices] = sort(bootstraps{iModel}.jointlikelihood,'descend');
        
        p = bootstraps{iModel};
        parameterEstimates.u0 = p.u0(:,mostLikelyIndices(1));
        parameterEstimates.v0 = p.v0(:,mostLikelyIndices(1));
        parameterEstimates.sigma_n = p.sigma_n(:,mostLikelyIndices(1));
        parameterEstimates.sigma_s = p.sigma_s(:,mostLikelyIndices(1));
        parameterEstimates.zeta = p.zeta(:,mostLikelyIndices(1));
        parameterEstimates.delta = p.delta(:,mostLikelyIndices(1));
        
        if shouldFitToSecondMomentOnly == 1
            q = x-mean(x,2);
            r = y-mean(y,2);
            [u_meso,v_meso,u_bg,v_bg,u_sm,v_sm,dmxdt,dmydt] = DecomposeTrajectories(q, r, t, parameterEstimates);
        else
            [u_meso,v_meso,u_bg,v_bg,u_sm,v_sm,dmxdt,dmydt] = DecomposeTrajectories(x, y, t, parameterEstimates);
        end
        
        % Confirm that the function being minimized does actually get
        % minimized when new parameters are added.
        u_residual = dmxdt - mean(u_meso,2);
        v_residual = dmydt - mean(v_meso,2);        
        min_sm = u_sm.^2 + v_sm.^2;
        min_cm = u_residual.^2 + v_residual.^2;
        epsilon2 = sum(sum(min_sm)) + sum(min_cm);
        u_meso_cm = u_meso - dmxdt;
        v_meso_cm = v_meso - dmydt;
        
        dqdt_meso = u_meso - mean(u_meso,2);
        drdt_meso = v_meso - mean(v_meso,2);
        
        FVUconst = FVU(u_sm,v_sm,u_sm+dqdt_meso,v_sm+drdt_meso);
        FDUconst = FDU(u_sm,v_sm,u_sm+dqdt_meso,v_sm+drdt_meso);
        

        cv_ms = dqdt_meso + sqrt(-1)*drdt_meso;
        cv_sm = u_sm + sqrt(-1)*v_sm;
        [psi,lambda]=sleptap(size(cv_ms,1));
        dt = t(2)-t(1);
        [f,sxx,syy,sxy]=mspec(dt,cv_ms,cv_sm,psi,'cyclic');
        gamma=frac(abs(sxy).^2,sxx.*syy);
        
        
        
        results{dof}.bootstraps{iModel} = bootstraps{iModel};
        results{dof}.bootstraps{iModel}.kappa = kappa(u_sm,v_sm);
        results{dof}.bootstraps{iModel}.FVU = FVUconst;
        results{dof}.bootstraps{iModel}.FDU = FDUconst;
        results{dof}.bootstraps{iModel}.epsilon2 = epsilon2;
        results{dof}.bootstraps{iModel}.epsilon2_sm = sum(sum(min_sm));
        results{dof}.bootstraps{iModel}.epsilon2_cm = sum(min_cm);
        results{dof}.bootstraps{iModel}.gamma = mean(mean(gamma(1:(floor(size(gamma,1)/2)),:)));
        
    end
    
end

fprintf('\n\\begin{tabular}{c|ccccccc}');
for dof = 1:6
    fprintf('\\multicolumn{8}{c}{}\\\\');
    fprintf('&  \\multicolumn{7}{c}{%d dof}\\\\',dof);
    fprintf('model & $\\kappa$ & FVU & FDU & $\\epsilon^2$ & $\\epsilon_{sm}^2$ & $\\epsilon_{cm}^2$ & $\\gamma$   \\\\ \\hline \n');
    bootstraps = results{dof}.bootstraps;
    for iModel = 1:length(bootstraps)
        p = bootstraps{iModel};
        fprintf('$%s$ & %.3f\t&%.3f\t&\t%.3f\t&\t%.3f&\t%.3f&\t%.3f&\t%.3f \\\\ \n',ModelParameter.modelName(p.parameters),p.kappa, p.FVU, p.FDU,p.epsilon2,p.epsilon2_sm,p.epsilon2_cm,p.gamma);
    end
end
fprintf('\\end{tabular}\n');

return
    
shouldShowFigure = 0;


load(sprintf('Rho%dDrifterSplineFits%d_K%d_dof%d.mat',SiteNumber,totalPermutations,K,dof));

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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% How many unique time series?
dt = 48; % window size
step = 1;
tStartIndex = (1:step:(length(t)-dt))';
tEndIndex = ((1+dt):step:length(t))';
tStartIndex = cat(1,tStartIndex,length(t)-dt);
tEndIndex = cat(1,tEndIndex,length(t));
nWindows = length(tStartIndex);
tFV = t(tStartIndex)+(t(tEndIndex)-t(tStartIndex))/2;

labels = cell(length(bootstraps),1);
targets = zeros(length(bootstraps),1);
if shouldShowFigure == 1
    figure('Position',[100 100 1000 800])
    top = subplot(3,1,1);
    middle = subplot(3,1,2);
    bottom = subplot(3,1,3);
end
for iModel = 1:length(bootstraps)

    [~,mostLikelyIndices] = sort(bootstraps{iModel}.jointlikelihood,'descend');
        
    p = bootstraps{iModel};
    u0 = p.u0(:,mostLikelyIndices(1));
    v0 = p.v0(:,mostLikelyIndices(1));
    sigma_n = p.sigma_n(:,mostLikelyIndices(1));
    sigma_s = p.sigma_s(:,mostLikelyIndices(1));
    zeta = p.zeta(:,mostLikelyIndices(1));
    delta = p.delta(:,mostLikelyIndices(1));
    
    u_meso = u0 + 0.5*(sigma_n+delta).*x + 0.5*(sigma_s-zeta).*y;
    v_meso = v0 + 0.5*(sigma_s+zeta).*x + 0.5*(delta-sigma_n).*y;
    
    u_bg = dmxdt - mean(u_meso,2);
    v_bg = dmydt - mean(v_meso,2);
    
    u_sm = dxdt - u_meso - u_bg;
    v_sm = dydt - v_meso - v_bg;
    
    FVUconst = FVU(u_sm,v_sm,dxdt-dmxdt,dydt-dmydt);
    FDUconst = FDU(u_sm,v_sm,dxdt-dmxdt,dydt-dmydt);
    FVUGlobalConst = FVU(dxdt-u_meso,dydt-v_meso,dxdt,dydt);
    
    fprintf('%s: %.3f\t&%.3f\t&\t%.3f\t&\t%.3f\n',ModelParameter.modelName(p.parameters),FVUGlobalConst,kappa(u_sm,v_sm) ,FVUconst,FDUconst  );
    if shouldShowFigure == 1
    FVUt = zeros(nWindows,1); FDUt = zeros(nWindows,1);  FVUGlobalt = zeros(nWindows,1);
    for iWindow = 1:nWindows
        tIndices = (tStartIndex(iWindow):tEndIndex(iWindow))';
        FVUt(iWindow) = FVU(u_sm(tIndices,:),v_sm(tIndices,:),dxdt(tIndices,:)-dmxdt(tIndices,:),dydt(tIndices,:)-dmydt(tIndices,:));
        FVUGlobalt(iWindow) = FVU(dxdt(tIndices,:)-u_meso(tIndices,:),dydt(tIndices,:)-v_meso(tIndices,:),dxdt(tIndices,:),dydt(tIndices,:));
        FDUt(iWindow) = FDU(u_sm(tIndices,:),v_sm(tIndices,:),dxdt(tIndices,:)-dmxdt(tIndices,:),dydt(tIndices,:)-dmydt(tIndices,:));
    end
    
    subplot(top)
    targets(iModel) = plot(tFV/86400,FVUGlobalt,'LineWidth',2); hold on
    h = gca;
    h.ColorOrderIndex= mod(h.ColorOrderIndex - 1 - 1,size(colororder,1))+1;
    plot([min(t) max(t)]/86400,FVUGlobalConst*[1 1],'LineWidth',1)
    
    subplot(middle)
    plot(tFV/86400,FVUt,'LineWidth',2); hold on
    h = gca;
    h.ColorOrderIndex= mod(h.ColorOrderIndex - 1 - 1,size(colororder,1))+1;
    plot([min(t) max(t)]/86400,FVUconst*[1 1],'LineWidth',1)
    
    subplot(bottom)
    plot(tFV/86400,FDUt,'LineWidth',2); hold on
    h = gca;
    h.ColorOrderIndex= mod(h.ColorOrderIndex - 1 - 1,size(colororder,1))+1;
    plot([min(t) max(t)]/86400,FDUconst*[1 1],'LineWidth',1)
    
    labels{iModel} = ModelParameter.modelName(p.parameters(2:end));
    end
end
if shouldShowFigure == 1
subplot(top)
legend(targets,labels)
xlim([min(t) max(t)]/86400)
ylabel('FVU Absolute');
title(sprintf('Site %d',SiteNumber));
legend(targets,labels)

subplot(middle)
xlim([min(t) max(t)]/86400)
ylabel('FVU');
subplot(bottom)
xlim([min(t) max(t)]/86400)
xlabel('t (days)');
ylabel('FDU');
packfig(3,1)
end


% print(sprintf('FVUFDURho%dDrifterSplineFits%d_K%d_dof%d.eps',SiteNumber,totalPermutations,K,dof),'-depsc2')