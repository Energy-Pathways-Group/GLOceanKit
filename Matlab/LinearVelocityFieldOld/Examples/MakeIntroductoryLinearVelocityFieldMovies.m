% 	Used to create movies showing single particle diffusivity
% 	(particles on a random walk) and how, with enough particles, that
% 	becomes a tracer diffusivity. I recommend a movie with 100
% 	particles, a movie with 500 particles, and then a movie with 500
% 	particles shown as a concentration.

scenario = 2;
showEndStateOnly = 1; % if 1, will show last time point, else will make a movie

% defaults
kappa = 0.5;
zeta = 0.0;
sigma = 0;
theta = 0;

if scenario == 1
    name = 'Diffusion';
	FramesFolder ='/Users/jearly/Desktop/diffusion_100_particles';
	kappa = 0.5;
	N=100;
	shouldShowTails = 0;
	shouldShowAsConcentration = 1;
    shouldShowVelocityField = 0;
elseif scenario == 2
    name = 'Diffusion';
	FramesFolder ='/Users/jearly/Desktop/diffusion_500_particles';
	kappa = 0.5;
	N=500;
	shouldShowTails = 0;
	shouldShowAsConcentration = 0;	
    shouldShowVelocityField = 0;
elseif scenario == 3
    name = 'Diffusion';
	FramesFolder ='/Users/jearly/Desktop/diffusion_concentration';
	kappa = 0.5;
	N=500;
	shouldShowTails = 0;
	shouldShowAsConcentration = 1;	
    shouldShowVelocityField = 0;
elseif scenario == 4
    name = 'Strain-Diffusion';
	FramesFolder ='/Users/jearly/Desktop/strain_diffusion_100_particles';
	kappa = 0.5;
	N=20;
    sigma = 4E-6;
    theta = -32*pi/180;
	shouldShowTails = 0;
	shouldShowAsConcentration = 0;
    shouldShowVelocityField = 1;
elseif scenario == 5
    name = 'Vorticity-Strain Dominated-Diffusion';
	FramesFolder ='/Users/jearly/Desktop/vorticity_strain_diffusion_100_particles';
	kappa = 0.5;
	N=20;
    s=2E-6;
    alpha = 1;
    sigma = 4E-6;
    zeta = 3.9E-6;
    theta = 0; %-32*pi/180;
	shouldShowTails = 0;
	shouldShowAsConcentration = 0;	
    shouldShowVelocityField = 1;
elseif scenario == 6
    name = 'Vorticity Dominated-Strain-Diffusion';
	FramesFolder ='/Users/jearly/Desktop/vorticity_dominated_strain_diffusion_100_particles';
	kappa = 0.5;
	N=20;
    s=2*2*pi/(7*86400);
    alpha = 1;
    sigma = abs(s*sinh(alpha));
    zeta = s*cosh(alpha);
    theta = -32*pi/180;
	shouldShowTails = 0;
	shouldShowAsConcentration = 0;	
    shouldShowVelocityField = 1;
elseif scenario == 7
    name = 'Vorticity-Diffusion';
	FramesFolder ='/Users/jearly/Desktop/vorticity_diffusion_100_particles';
	kappa = 0.5;
	N=20;
    sigma = 0;
    zeta = 4E-6;
    theta = 0;
	shouldShowTails = 0;
	shouldShowAsConcentration = 0;
    shouldShowVelocityField = 1;
elseif scenario == 8
    name = 'Vorticity-Strain-Matched-Diffusion';
	FramesFolder ='/Users/jearly/Desktop/vorticity_strain_matched_diffusion_100_particles';
	kappa = 0.5;
	N=20;
    sigma = 4.0E-6;
    zeta = 4E-6;
    theta = 0;
	shouldShowTails = 0;
	shouldShowAsConcentration = 0;
    shouldShowVelocityField = 1;
end

shouldShowTheoreticalSecondMoment = 1;

t=(0:30*60:7*86400)';
sigma_n = sigma*cos(2*theta);
sigma_s = sigma*sin(2*theta);
params = LinearVelocityField(zeta,sigma,theta,kappa);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Create the initial float positions
%

x=linspace(-1000,1000,N);
y=linspace(-1000,1000,N);
[X,Y]=meshgrid(x,y);
x0 = reshape(X,1,N*N);
y0 = reshape(Y,1,N*N);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Create the time series
%

[x, y] = ParticlePathInLinearVelocityField( x0, y0, t, params );


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% 	Make the frames folder
%
if showEndStateOnly == 0
    if exist(FramesFolder, 'dir') == 0
        mkdir(FramesFolder);
    end
end

numDrifters = size(x,2);
days = t/86400;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Center of mass coordinates
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[x_com, y_com, q, r] = CenterOfMass( x, y );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Axis range
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nbins=50;
maxrange = 1.1*max(max(abs([q;r])));
xrange = [-maxrange maxrange];
yrange = [-maxrange maxrange];

Area_per_bin = (max(xrange)-min(xrange))*(max(yrange)-min(yrange))/(nbins*nbins);
Area_per_drifter = (max(x0)-min(x0))*(max(y0)-min(y0))/(N*N);
drifters_per_bin = Area_per_bin / Area_per_drifter;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Strain model velocity field
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

npoints = 10;
aspect_ratio = (xrange(end)-xrange(1))/(yrange(end)-yrange(1));
xdim = linspace(xrange(1), xrange(end), npoints);
ydim = linspace(yrange(1), yrange(end), npoints/aspect_ratio);
[Y,X]=meshgrid(ydim,xdim);
X = reshape(X, 1, numel(X));
Y = reshape(Y, 1, numel(Y));
[u, v] = VelocityFromPositionInStrainVorticityField( X, Y, zeta, sigma_n, sigma_s );

X = reshape(X, length(xdim), length(ydim));
Y = reshape(Y, length(xdim), length(ydim));
u = reshape(u, length(xdim), length(ydim));
v = reshape(v, length(xdim), length(ydim));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Second moment matrix (from observations)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if shouldShowTheoreticalSecondMoment == 1
    [M_qq, M_rr, M_qr] = SecondMomentMatrix( q, r );
    [M_qq, M_rr, M_qr] = MomentTensorEvolutionInStrainVorticityField( M_qq(1), M_rr(1), M_qr(1), t, zeta, sigma, theta, kappa );
    theta_theoretical = zeros(size(t));
    minD_theoretical = zeros(size(t));
    maxD_theoretical = zeros(size(t));
    for iTime =1:length(t)
        [A, B] = eig([M_qq(iTime), M_qr(iTime); M_qr(iTime), M_rr(iTime)]);
        v_theo = A(:,2);
        theta_theoretical(iTime) = atan(v_theo(2)/v_theo(1));
        minD_theoretical(iTime) = B(1,1);
        maxD_theoretical(iTime) = B(end,end);
    end
    [k_theoretical,l_theoretical] = ab2kl(sqrt(maxD_theoretical),sqrt(minD_theoretical));
end

[minD, maxD, theta] = SecondMomentMatrix( q, r, 'eigen' );

% The extra sqrt(2) factor is the correct scaling if we wish to draw the e-folding scale
% of a Gaussian distribution. But, we don't!
[k,l] = ab2kl(sqrt(maxD),sqrt(minD));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Draw the drifter paths
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure('Position', [50 50 680 680])
set(gcf,'PaperPositionMode','auto')
set(gcf, 'Color', 'w');

rgb0=[1.0 1.0 1.0];
rgbN=[1.0 0.0 0.0];
rgb = interp1(0:1,[rgb0; rgbN],0:0.01:1.0);
colormap(rgb);

if showEndStateOnly == 1
    iTime0 = length(t);
else
    iTime0 = 1;
end

for iTime=iTime0:length(t)
	clf
    
    if (shouldShowVelocityField == 1)
        quiver(X, Y, u, v, 'LineWidth', 0.5, 'Color', 'b')
    end
    
	if (shouldShowTails == 1)
		% Determine how much of the path we want to plot
		tailLength=1440/4;
		tailLength=20;
		startIndex = iTime-tailLength;
		if (startIndex < 1) startIndex=1; end
		range=startIndex:iTime;
	

		if iTime > 1
			plot( q(range,:), r(range,:), 'LineWidth', 2)
		end
	end
	
	hold on
	
    if (shouldShowAsConcentration == 1)
        [c,x,y] = twodhist( q(iTime,:)', r(iTime,:)', xrange, yrange, nbins );
        c = c / drifters_per_bin;
        
        pcolor(x-(x(2)-x(1))/2,y-(y(2)-y(1))/2,c); shading flat; caxis([0 1.0]);
    else
        scatter( q(iTime,:), r(iTime,:), 8^2, 'r', 'fill')
    end
    
    if shouldShowTheoreticalSecondMoment == 1
        ellipseplot(k_theoretical(iTime),l_theoretical(iTime),theta_theoretical(iTime),[0,0], '5b-')
    end
    
	ellipseplot(k(iTime),l(iTime),theta(iTime),[0,0], '3k-')
    

	
 	set(gca,'FontSize', 16)
 	if (shouldShowAsConcentration == 1)
		title( sprintf('Drifter concentration on day %d at %2d:%02d hours', floor(t(iTime)/86400), floor(mod(t(iTime)/3600,24)), floor(mod(t(iTime)/60,60)) ), 'fontsize', 24, 'FontName', 'Helvetica' );
 	else
 	 	title( sprintf('Drifter positions on day %d at %2d:%02d hours', floor(t(iTime)/86400), floor(mod(t(iTime)/3600,24)), floor(mod(t(iTime)/60,60)) ), 'fontsize', 24, 'FontName', 'Helvetica' );
 	end
	xlabel('meters', 'FontSize', 20.0, 'FontName', 'Helvetica');
	ylabel('meters', 'FontSize', 20.0, 'FontName', 'Helvetica');
	axis([xrange(1) xrange(2)  yrange(1) yrange(2)])
    axis equal
		
    if showEndStateOnly == 1
        continue;
    end

	% write everything out
	output = sprintf('%s/t_%03d', FramesFolder,iTime-1);
	ScaleFactor=4;
	if (shouldShowAsConcentration == 1)
		print(sprintf('-r%d',72*ScaleFactor), '-dpng', output );
	else
		print('-depsc2', output)
	end
end