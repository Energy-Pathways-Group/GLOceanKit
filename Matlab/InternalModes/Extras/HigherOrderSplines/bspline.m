%bspline
%   
function X = bspline( t, t_knot, K )

dt_knot = diff(t_knot);
if ~isempty(find(dt_knot <= 0,1))
    disp('t_knot must be strictly monotonically increasing');
    return
end


S = K-1;

% Now we increase the multiplicity of the knot points at the beginning and
% the end of the interval so that the splines do not extend past the end
% points.
t_knot = [repmat(t_knot(1),S,1); t_knot; repmat(t_knot(end),S,1)];
dt_knot = diff(t_knot);
t_knot2 = t_knot + [dt_knot; 0];

% numer of knots
M = length(t_knot);

% This is true assuming the original t_knot was strictly monotonically
% increasing (no repeat knots) and we added repeat knots at the beginning
% and end of the sequences.
N_splines = M - K;

% number of collocation points
N = length(t);

% Rows are the N collocation points
% Columns are the M splines
X = zeros(N,N_splines,K); % This will contain all splines and their derivatives
XB = zeros(N,N_splines,K); % This will contain all splines through order K
for t_i=1:N % loop through all N collocation points
    i = find( t_knot <= t(t_i) & t(t_i) < t_knot2, 1, 'last' );
    if isempty(i)
        if t(t_i) < t_knot(1)
%             i = find( t_knot == t_knot(1), 1, 'last');
            continue; %This continue means we don't need to set b(1) = 0; or check indices on the delta_r line
        elseif t(t_i) == t_knot(end)
            i = find( t_knot < t(t_i), 1, 'last'); 
        else
%             i = find( t_knot < t_knot(end), 1, 'last');
            continue; %b(1) = 0;
        end
    end
    
    delta_r = zeros(K,1);
    delta_l = zeros(K,1);
    
    XB(t_i,i,1) = 1;
    
    b = zeros(K,1); b(1) = 1;
    for j=1:(K-1) % loop through splines of increasing order
       delta_r(j) = t_knot(i+j) - t(t_i);
       delta_l(j) = t(t_i) - t_knot(i+1-j);
       
       saved = 0;
       for r=1:j % loop through the nonzero splines
           term = b(r)/(delta_r(r) + delta_l(j+1-r));
           b(r) = saved + delta_r(r)*term;
           saved = delta_l(j+1-r)*term;
       end
       b(j+1) = saved;
       
       indices = max(1,i-j):i;
       XB(t_i,indices,j+1) = b(1:length(indices));
    end
    
    indices = max(1,i-K+1):i;
    X(t_i,indices,1) = b(1:length(indices));
    
end

diff_coeff = @(a,r,m) (K-m)*(a(2)-a(1))/(t_knot(r+K-m) - t_knot(r));

for r=1:N_splines
    % alpha mimics equation X.16 in deBoor's PGS, but localized to avoid
    % the zero elements.
    alpha = zeros(S+2,S+2); % row is the coefficient, column is the derivative (1=0 derivatives)
    alpha(2,1) = 1;
    for m=1:S
        for i=1:(m+1)
            a = alpha(:,m);
            alpha(i+1,m+1) = diff_coeff(a(i:end),r+i-1,m);
            if isinf(alpha(i+1,m+1)) || isnan(alpha(i+1,m+1))
                alpha(i+1,m+1) = 0;
            end
            if r+i-1>N_splines
                B = zeros(N,1);
            else
                B = XB(:,r+i-1,K-m);
            end
            X(:,r,m+1) = X(:,r,m+1) + alpha(i+1,m+1)*B;
        end
    end
end
