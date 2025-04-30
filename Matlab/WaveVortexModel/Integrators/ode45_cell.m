function [T,Y] = ode45_cell(odefun,tspan,y0,options)
%ODE45_CELL  Dormand–Prince (4,5) ODE solver with *cell-array* states.
%
%   [T,Y] = ODE45_CELL(ODEFUN,TSPAN,Y0,OPTIONS) behaves like MATLAB®
%   ODE45 but accepts Y0 and returns Y as *cell arrays*.  Each cell may
%   contain any numeric array (scalar → N‑D).
%
%   • TSPAN:  [t0 tf]  or a vector of desired output times.  If a vector is
%             supplied the solver lands exactly on those times and calls
%             OutputFcn there.
%   • OPTIONS fields (all optional)
%         RelTol            — scalar (default 1e‑3)
%         AbsTol | absTolerance — scalar or cell (default 1e‑6)
%         InitialStep       — scalar first step • MaxStep — scalar limit
%         OutputFcn(t,y,flag)   where flag = 'init','', or 'done'
%
%   All other semantics follow MATLAB ODE45.
% -------------------------------------------------------------------------

% ─── Dormand–Prince tableau ─────────────────────────────────────────────
A = [0 0 0 0 0 0;
     1/5 0 0 0 0 0;
     3/40 9/40 0 0 0 0;
     44/45 -56/15 32/9 0 0 0;
     19372/6561 -25360/2187 64448/6561 -212/729 0 0;
     9017/3168 -355/33 46732/5247 49/176 -5103/18656 0];
B5 = [35/384 0 500/1113 125/192 -2187/6784 11/84 0];   % 5th order
B4 = [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40];
C  = [0; 1/5; 3/10; 4/5; 8/9; 1];
order = 5;

%% ─── Parse inputs ──────────────────────────────────────────────────────
if nargin<4, options = struct(); end
if ~iscell(y0), error('y0 must be a cell array'); end
if ~isa(odefun,'function_handle'), error('odefun must be a function handle'); end

isVectorTspan = numel(tspan) > 2;
if ~isVectorTspan && numel(tspan)~=2
    error('tspan must be [t0 tf] or a monotone vector'); end

RelTol      = getOpt(options,'RelTol',1e-3);
AbsTolIn    = getOpt(options,'AbsTol',[]);
AbsTolInAlt = getOpt(options,'absTolerance',[]);
if isempty(AbsTolIn), AbsTolIn = AbsTolInAlt; end
if isempty(AbsTolIn), AbsTolIn = 1e-6; end
InitialStep = getOpt(options,'InitialStep',[]);
MaxStep     = getOpt(options,'MaxStep',max(tspan)-min(tspan));
OutputFcn   = getOpt(options,'OutputFcn',[]);

% Expand AbsTol to cell
if ~iscell(AbsTolIn)
    AbsTol = cellfun(@(y) repmat(AbsTolIn,size(y)), y0,'UniformOutput',false);
else
    if ~isequal(size(AbsTolIn),size(y0))
        error('AbsTol/absTolerance cell array must match y0'); end
    AbsTol = AbsTolIn;
end

%% ─── Allocate storage ──────────────────────────────────────────────────
if isVectorTspan
    T = tspan(:);
    nOut = numel(T);
    Y = cellfun(@(y) zeros([nOut,size(y)]), y0,'UniformOutput',false);
    for k=1:numel(Y), Y{k}(1,:) = y0{k}(:).'; end
    nextOut = 2;
else
    alloc = 1000;
    T = zeros(alloc,1); T(1)=tspan(1);
    Y = cellfun(@(y) zeros([alloc,size(y)]), y0,'UniformOutput',false);
    for k=1:numel(Y), Y{k}(1,:) = y0{k}(:).'; end
    idxStore = 1;
end

%% ─── Initial step size ────────────────────────────────────────────────
h = initialStep(odefun,tspan(1),y0,RelTol,AbsTol,InitialStep);

if ~isempty(OutputFcn)
    if feval(OutputFcn,tspan(1),y0,'init'), return; end
end

t = tspan(1);

%% ─── Main loop ────────────────────────────────────────────────────────
while true
    if isVectorTspan
        if nextOut>numel(T), break; end
        h = min([h, MaxStep, T(nextOut)-t]);
    else
        h = min([h, MaxStep, tspan(2)-t]);
    end
    if h <= eps(max(1,abs(t)))
        warning('ode45_cell:stepUnderflow','Step size underflow.'); break; end

    % --- Stage evaluations (k1..k6) -----------------------------------
    K1 = odefun(t,                    y0);
    K2 = odefun(t+C(2)*h, add( y0, h*A(2,1), K1));
    K3 = odefun(t+C(3)*h, add3(y0, h*A(3,1), K1, h*A(3,2), K2));
    K4 = odefun(t+C(4)*h, add4(y0, h*A(4,1), K1, h*A(4,2), K2, h*A(4,3), K3));
    K5 = odefun(t+C(5)*h, add5(y0, h*A(5,1), K1, h*A(5,2), K2, h*A(5,3), K3, h*A(5,4), K4));
    K6 = odefun(t+C(6)*h, add6(y0, h*A(6,1), K1, h*A(6,2), K2, h*A(6,3), K3, h*A(6,4), K4, h*A(6,5), K5));

    y5 = add7(y0, h*B5(1),K1, h*B5(3),K3, h*B5(4),K4, h*B5(5),K5, h*B5(6),K6);
    y4 = add7(y0, h*B4(1),K1, h*B4(3),K3, h*B4(4),K4, h*B4(5),K5, h*B4(6),K6, h*B4(7),K6);

    err = errorNorm(y4,y5,AbsTol,RelTol);

    if err <= 1      % Accept step ------------------------------------
        t  = t + h;
        y0 = y5;
        if isVectorTspan
            if abs(t - T(nextOut)) <= 1e-12*max(1,abs(t))
                for k=1:numel(Y), Y{k}(nextOut,:) = y0{k}(:).'; end
                if ~isempty(OutputFcn)
                    if feval(OutputFcn,t,y0,''), break; end
                end
                nextOut = nextOut + 1;
            end
            if t >= T(end)-eps(max(1,abs(T(end)))), break; end
        else
            idxStore = idxStore + 1;
            if idxStore>numel(T)
                T = [T; zeros(1000,1)];
                Y = cellfun(@(Yi) cat(1,Yi,zeros([1000,size(Yi,2:size(Yi,2))])), Y,'UniformOutput',false);
            end
            T(idxStore)=t;
            for k=1:numel(Y), Y{k}(idxStore,:) = y0{k}(:).'; end
            if ~isempty(OutputFcn)
                if feval(OutputFcn,t,y0,''), break; end
            end
            if t >= tspan(2)-eps(max(1,abs(tspan(2)))), break; end
        end
    end

    % --- Step-size update ---------------------------------------------
    h = 0.9*h*max(0.2,min(5, err^(-1/(order+1))));
end

if ~isempty(OutputFcn), feval(OutputFcn,t,y0,'done'); end

if ~isVectorTspan
    T = T(1:idxStore);
    Y = cellfun(@(Yi) Yi(1:idxStore,:), Y,'UniformOutput',false);
end

end  % ode45_cell

%% Helper utilities -------------------------------------------------------
function val=getOpt(s,f,d); if isfield(s,f), val=s.(f); else, val=d; end; end

% Vectorised cell additions ---------------------------------------------
function y = add(y0,a1,k1)
    y = cellfun(@(y_,k_) y_ + a1*k_, y0,k1,'UniformOutput',false); end
function y = add3(y0,a1,k1,a2,k2)
    y = cellfun(@(y_,k1_,k2_) y_ + a1*k1_ + a2*k2_, y0,k1,k2,'UniformOutput',false); end
function y = add4(y0,a1,k1,a2,k2,a3,k3)
    y = cellfun(@(y_,k1_,k2_,k3_) y_ + a1*k1_ + a2*k2_ + a3*k3_, y0,k1,k2,k3,'UniformOutput',false); end
function y = add5(y0,a1,k1,a2,k2,a3,k3,a4,k4)
    y = cellfun(@(y_,k1_,k2_,k3_,k4_) y_ + a1*k1_ + a2*k2_ + a3*k3_ + a4*k4_, y0,k1,k2,k3,k4,'UniformOutput',false); end
function y = add6(y0,a1,k1,a2,k2,a3,k3,a4,k4,a5,k5)
    y = cellfun(@(y_,k1_,k2_,k3_,k4_,k5_) y_ + a1*k1_ + a2*k2_ + a3*k3_ + a4*k4_ + a5*k5_, y0,k1,k2,k3,k4,k5,'UniformOutput',false); end
function y = add7(y0,a1,k1,a3,k3,a4,k4,a5,k5,a6,k6,a7,k7)
    if nargin<15, a7=0; k7=k6; end
    y = cellfun(@(y_,k1_,k3_,k4_,k5_,k6_,k7_) y_ + a1*k1_ + a3*k3_ + a4*k4_ + a5*k5_ + a6*k6_ + a7*k7_, y0,k1,k3,k4,k5,k6,k7,'UniformOutput',false); end

% Error norm --------------------------------------------------------------
function err = errorNorm(y4,y5,AbsTol,RelTol)
    denom = cellfun(@(s4,s5,a) a + max(abs(s4),abs(s5))*RelTol, y4,y5,AbsTol,'UniformOutput',false);
    E     = cellfun(@(s4,s5,d) (s5-s4)./d, y4,y5,denom,'UniformOutput',false);
    err   = sqrt(sum(cellfun(@(e) sum(e(:).^2), E))/numel(E{1}(:))); end

% Initial step heuristic --------------------------------------------------
function h = initialStep(odefun,t0,y0,RelTol,AbsTol,hUser)
    if ~isempty(hUser), h = hUser; return; end
    f0 = odefun(t0,y0);
    scale = cellfun(@(y,a) a + abs(y)*RelTol, y0,AbsTol,'UniformOutput',false);
    d0 = sqrt(sum(cellfun(@(f,s) sum((f(:)./s(:)).^2), f0,scale))/numel(scale{1}(:)));
    h = 0.01*max(1e-6, 1/d0^(1/(order+1))); end
