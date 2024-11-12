N = 10;
T = 1;
t=T*(0:(N-1))'/(N-1);
dt = t(2)-t(1);
df = 1/(2*(N-1)*dt);

t2 = linspace(t(1),t(end),1000).';

D = CosineTransformBackMatrix(N);

%% This is just the requirement that \int F dz = known_value
b = zeros(N,1);
b(1) = 1;
w = D.'\b;

x = D(:,N-1);
(1/T)*(sum(x(2:end-1).*x(2:end-1))*dt+x(1)*x(1)*dt/2 + x(end)*x(end)*dt/2)

norm = 0.5*ones(N,1);
norm(1) = 0.25;
norm(end) = 1.0;

% This does not work. Singular matrix, suggesting it is not sufficiently
% constrained
Q = (D.').^2;
% w = Q \ norm;

w_actual = ones(N,1);
w_actual(1) = 0.5;
w_actual(end) = 0.5;
w_actual = w_actual/(N-1);

Q*w_actual

%%
Nz = 10;
Nj = Nz-1;
Lz = 1300;
N0 = 5.2e-3;
N2 = @(z) N0*N0*exp(2*z/L_gm);
N2 = @(z) N0*N0*ones(size(z));
z = WVStratifiedFlow.quadraturePointsForStratifiedFlow(Lz,Nz,N2=N2,latitude=33);
verticalModes = InternalModesSpectral(N2=N2,zIn=[-Lz 0],zOut=z,latitude=33,nModes=Nj,nEVP=max(256,floor(2.1*Nz)));
verticalModes.normalization = Normalization.geostrophic;
verticalModes.upperBoundary = UpperBoundary.rigidLid;
[Finv,Ginv,h,k,int_n2g,int_f,int_g] = verticalModes.ModesAtFrequency(0,'int_N2_G_dz/g','int_F_dz','int_G_dz');
[P,Q,PFinv,PF,QGinv,QG,h,w] = WVStratifiedFlow.verticalProjectionOperatorsWithRigidLid(Finv,Ginv,h,Nj,Lz);
int_n2g=int_n2g.'; int_f=int_f.'; int_g=int_g.';

% probably don't want to chop this mode off
int_f = cat(1,Lz,int_f(1:end-1));
w = PFinv.'\int_f

% PFinv * (P.*A0)

%%
w = zeros(size(t));
for i=1:length(w)
    Ci = cardinalFunctionSine(t,i,N);
    % Ci = cardinalFunction(t,i);
    w(i) = integral( Ci, t(1), t(end));
end

function Ci = cardinalFunctionSine(xj,i,N)
xi = xj(i);
Ci = @(x) sin((N-1)*pi*x)./((N-1)*pi*cos((N-1)*pi*xi).*(x-xi));
end

function Ci = cardinalFunction(xj,i)
xi = xj(i);
xj(i) = [];
denom = prod(xi - xj);
Ci = @(x) prod(x-xj)/denom;
end

% w = Q \ norm;