SiteNumber=1;
load(sprintf('smoothedGriddedRho%dDrifters.mat',SiteNumber));
% The last drifter at both sites is only partial time series.
x = x(:,1:(end-1));
y = y(:,1:(end-1));
nDrifters = size(x,2);
nT = size(x,1);

shouldFitToSecondMomentOnly = 0;

FramesFolder = './BootstrapData';
if exist(FramesFolder,'dir') == 0
	mkdir(FramesFolder);
end

% How many bootstrap combinations?
totalPermutations = 1000;

for dof = 1:6
    
    if shouldFitToSecondMomentOnly == 1
        filename = sprintf('BootstrapData/Rho%dDrifterSpline2ndMomFits%d_dof%d.mat',SiteNumber,totalPermutations,dof);
        u0v0 = [];
    else
        filename = sprintf('BootstrapData/Rho%dDrifterSplineFits%d_dof%d.mat',SiteNumber,totalPermutations,dof);
        u0v0 = ModelParameter.u0v0;
    end
    % Pre-allocate our arrays
    drifterIDs = cell(totalPermutations,1);
    nothing = zeros(nT,totalPermutations);
    emptyStruct = struct('parameters',[],'dof',dof,'u0', nothing,'v0', nothing,'sigma_n', nothing, 'sigma_s', nothing, 'zeta', nothing,'delta', nothing);
    
    b = 0;
    clear bootstraps;
    % b=b+1; bootstraps{b} = emptyStruct; bootstraps{b} = [ModelParameter.u0v0];
    b=b+1; bootstraps{b} = emptyStruct; bootstraps{b}.parameters = [u0v0,ModelParameter.vorticity];
    b=b+1; bootstraps{b} = emptyStruct; bootstraps{b}.parameters = [u0v0,ModelParameter.divergence];
    b=b+1; bootstraps{b} = emptyStruct; bootstraps{b}.parameters = [u0v0,ModelParameter.strain];
    b=b+1; bootstraps{b} = emptyStruct; bootstraps{b}.parameters = [u0v0,ModelParameter.vorticity,ModelParameter.divergence];
    b=b+1; bootstraps{b} = emptyStruct; bootstraps{b}.parameters = [u0v0,ModelParameter.strain,ModelParameter.vorticity];
    b=b+1; bootstraps{b} = emptyStruct; bootstraps{b}.parameters = [u0v0,ModelParameter.strain,ModelParameter.divergence];
    b=b+1; bootstraps{b} = emptyStruct; bootstraps{b}.parameters = [u0v0,ModelParameter.strain,ModelParameter.vorticity,ModelParameter.divergence];
    if dof == 1 && shouldFitToSecondMomentOnly == 0
        b=b+1; bootstraps{b} = emptyStruct; bootstraps{b}.parameters = [u0v0,ModelParameter.u1v1,ModelParameter.vorticity];
        b=b+1; bootstraps{b} = emptyStruct; bootstraps{b}.parameters = [u0v0,ModelParameter.u1v1,ModelParameter.divergence];
        b=b+1; bootstraps{b} = emptyStruct; bootstraps{b}.parameters = [u0v0,ModelParameter.u1v1,ModelParameter.strain];
        b=b+1; bootstraps{b} = emptyStruct; bootstraps{b}.parameters = [u0v0,ModelParameter.u1v1,ModelParameter.vorticity,ModelParameter.divergence];
        b=b+1; bootstraps{b} = emptyStruct; bootstraps{b}.parameters = [u0v0,ModelParameter.u1v1,ModelParameter.strain,ModelParameter.vorticity];
        b=b+1; bootstraps{b} = emptyStruct; bootstraps{b}.parameters = [u0v0,ModelParameter.u1v1,ModelParameter.strain,ModelParameter.divergence];
        b=b+1; bootstraps{b} = emptyStruct; bootstraps{b}.parameters = [u0v0,ModelParameter.u1v1,ModelParameter.strain,ModelParameter.vorticity,ModelParameter.divergence];
    end
    
    % Now let's walk through all the possible combinations
    drifterPermutations = zeros(totalPermutations,nDrifters);
    
    tic
    for iPermutation = 1:totalPermutations
        drifterIndices = 1+floor(rand(1,nDrifters)*nDrifters);
        drifterPermutations(iPermutation,:) = drifterIndices;
        if mod(iPermutation,50) == 0
            fprintf('%d of %d permutations.\n',iPermutation,totalPermutations);
        end
        
        t_sub = t;
        x_sub = x(:,drifterIndices);
        y_sub = y(:,drifterIndices);
        
        if shouldFitToSecondMomentOnly == 1
            x_sub = x_sub - mean(x_sub,2);
            y_sub = y_sub - mean(y_sub,2);
        end
        
        for iModel = 1:length(bootstraps)
            [p,B] = EstimateLinearVelocityFieldParameters( x_sub, y_sub, t_sub, bootstraps{iModel}.parameters, bootstraps{iModel}.dof);
            if dof == 1
                bootstraps{iModel}.u0(:,iPermutation)=B*p.u0;
                bootstraps{iModel}.v0(:,iPermutation)=B*p.v0;
                bootstraps{iModel}.u1(:,iPermutation)=B*p.u1;
                bootstraps{iModel}.v1(:,iPermutation)=B*p.v1;
                bootstraps{iModel}.sigma_n(:,iPermutation)=B*p.sigma_n;
                bootstraps{iModel}.sigma_s(:,iPermutation)=B*p.sigma_s;
                bootstraps{iModel}.zeta(:,iPermutation)=B*p.zeta;
                bootstraps{iModel}.delta(:,iPermutation)=B*p.delta;
            else
                bootstraps{iModel}.u0(:,iPermutation)=p.u0;
                bootstraps{iModel}.v0(:,iPermutation)=p.v0;
                bootstraps{iModel}.u1(:,iPermutation)=p.u1;
                bootstraps{iModel}.v1(:,iPermutation)=p.v1;
                bootstraps{iModel}.sigma_n(:,iPermutation)=p.sigma_n;
                bootstraps{iModel}.sigma_s(:,iPermutation)=p.sigma_s;
                bootstraps{iModel}.zeta(:,iPermutation)=p.zeta;
                bootstraps{iModel}.delta(:,iPermutation)=p.delta;
            end
        end
    end
    toc
    
    dt = 6;
    for iModel = 1:length(bootstraps)
        bootstraps{iModel}.jointlikelihood = EstimateSolutionLikelihoodFromBootstraps(bootstraps{iModel},bootstraps{iModel}.parameters,dt);
    end
    
    save(filename,'bootstraps','drifterPermutations','B');
end

return;

m = bootstraps_strain_vorticity;
jointlikelihood = EstimateSolutionLikelihoodFromBootstraps(m,m.parameters,dt);
[val,idx] = max(jointlikelihood);
[~,mostLikelyIndices] = sort(jointlikelihood,'descend');


figure
subplot(2,1,1)
plot(t/86400,m.sigma(:,mostLikelyIndices(1:900)),'LineWidth',1,'Color',0.5*[1 0 0]), hold on
ylabel('\sigma')
plot(t/86400,m.sigma(:,idx),'LineWidth',2,'Color',0*[1 1 1])

subplot(2,1,2)
plot(t/86400,m.theta(:,mostLikelyIndices(1:900))*180/pi,'LineWidth',1,'Color',0.5*[1 0 0]), hold on
plot(t/86400,m.theta(:,idx)*180/pi,'LineWidth',2,'Color',0*[1 1 1])
ylabel('\theta')
% print('SplineFitSketch.png','-dpng');

return 
m = b;
% m = cm_model_strain;
iTime = 5*48;
[bandwidth,density,X,Y]=kde2d(cat(2,m.sigma_n(iTime,:).',m.sigma_s(iTime,:).'));
pctTarget = flip(0.1:0.1:0.8);
    
    dLevels = DensityLevelForCDF(X,Y,density, pctTarget);
    
    pctLabels = cell(length(pctTarget),1);
    for i=1:length(pctTarget)
        pctLabels{i} = sprintf('%d %%',round(pctTarget(i)*100));
    end
    
figure
M = contourf(X,Y,density,dLevels);
hold on, scatter(m.sigma_n(iTime,:),m.sigma_s(iTime,:),2^2,.5*[1 1 1])
hold on
for i=2:2:10
    rectangle('Position',[-i -i 2*i 2*i]*1e-6, 'Curvature', [1 1]);
end
axis equal
    xlabel('\sigma_n')
    ylabel('\sigma_s')
cmap = colormap;
cmap(1,:)=1;
colormap(cmap)
cb = colorbar('eastoutside');
    cb.Ticks = dLevels;
    cb.TickLabels = pctLabels;
xlim([-1e-5 1e-5])
ylim([-1e-5 1e-5])