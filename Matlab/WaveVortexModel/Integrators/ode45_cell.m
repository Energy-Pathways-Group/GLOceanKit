function [T,Y] = ode45_cell(odefun,tspan,y0,options)
%ODE45_CELL  Dormand–Prince (4,5) ODE solver that works with *cell arrays*.
%
%   This routine mirrors the behaviour of MATLAB® ODE45 but allows the
%   initial condition Y0, the running solution Y, and the absolute
%   tolerances to be given as *cell arrays*.  Each cell may contain a
%   scalar, vector, matrix, or N‑D numeric array – including complex
%   values.
%
%   Syntax
%   ------
%       [T,Y] = ode45_cell(ODEFUN,TSPAN,Y0)
%       [T,Y] = ode45_cell(ODEFUN,TSPAN,Y0,OPTIONS)
%
%   Inputs
%   ------
%   ODEFUN   – function handle @(t,y) that returns dy/dt in the same cell
%               layout as Y0.
%   TSPAN    – Either a 2‑element vector [t0 tf] *or* a monotonically
%               increasing (or decreasing) vector of output times.  If a
%               full vector is supplied, the solver stops exactly at those
%               points.
%   Y0       – Cell array of initial state values.
%   OPTIONS  – Struct with optional fields (default shown in brackets):
%                 RelTol          (1e‑3)
%                 AbsTol | absTolerance   (1e‑6)
%                 InitialStep     ([]   auto)
%                 MaxStep         (tf‑t0)
%                 OutputFcn@(t,y,state)  where state ∈ {'init','', 'done'}
%
%   Outputs
%   -------
%   T   – Column vector of times.
%   Y   – Cell array, size(Y) == size(Y0).  Y{i} has leading dimension
%         length(T) and trailing dims equal to size(Y0{i}).
%
%   ---------------------------------------------------------------------

%% Dormand–Prince (4,5) Butcher tableau ---------------------------------
A = [0 0 0 0 0 0;
     1/5 0 0 0 0 0;
     3/40 9/40 0 0 0 0;
     44/45 -56/15 32/9 0 0 0;
     19372/6561 -25360/2187 64448/6561 -212/729 0 0;
     9017/3168 -355/33 46732/5247 49/176 -5103/18656 0];
B5 = [35/384 0 500/1113 125/192 -2187/6784 11/84 0];      % fifth‑order
B4 = [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40]; % fourth‑order
C  = [0; 1/5; 3/10; 4/5; 8/9; 1];
METHOD_ORDER = 5;

%% ---------------- Input checking & defaults ---------------------------
if nargin < 4, options = struct(); end
if ~iscell(y0)
    error('Y0 must be a cell array.');
end
if ~isa(odefun,'function_handle')
    error('ODEFUN must be a function handle.');
end
% tspan validity
isVectorTspan = numel(tspan) > 2;
if ~isVectorTspan && numel(tspan) ~= 2
    error('TSPAN must be a 2‑element interval or a monotone vector.');
end

RelTol      = getOpt(options,'RelTol',1e-3);
AbsTolIn    = getOpt(options,'AbsTol',[]);
if isempty(AbsTolIn), AbsTolIn = getOpt(options,'absTolerance',1e-6); end
if isempty(AbsTolIn), AbsTolIn = 1e-6; end
InitialStep = getOpt(options,'InitialStep',[]);
MaxStep     = getOpt(options,'MaxStep',max(tspan)-min(tspan));
OutputFcn   = getOpt(options,'OutputFcn',[]);

% Expand AbsTol to a cell array
if ~iscell(AbsTolIn)
    AbsTol = cellfun(@(y) repmat(AbsTolIn,size(y)), y0,'UniformOutput',false);
else
    if ~isequal(size(AbsTolIn),size(y0))
        error('AbsTol cell size must match Y0.');
    end
    AbsTol = AbsTolIn;
end

%% -------------------- Pre‑allocate output -----------------------------
if isVectorTspan
    T = tspan(:);
    nOut = numel(T);
    Y = cellfun(@(y) zeros([nOut, size(y)]), y0,'UniformOutput',false);
    for i = 1:numel(Y), Y{i}(1,:) = y0{i}(:).'; end
    nextOut = 2;                           % index of next requested time
else
    chunk = 1000;
    T = zeros(chunk,1); T(1) = tspan(1);
    Y = cellfun(@(y) zeros([chunk, size(y)]), y0,'UniformOutput',false);
    for i = 1:numel(Y), Y{i}(1,:) = y0{i}(:).'; end
    storeIdx = 1;                          % last filled index
end

%% -------------------- Initial step size ------------------------------
h = initialStep(odefun,tspan(1),y0,RelTol,AbsTol,InitialStep,METHOD_ORDER);

t = tspan(1);

% Call OutputFcn with ''init''
if ~isempty(OutputFcn)
    if feval(OutputFcn,t,y0,'init'), return; end
end

%% ==================== Main integration loop ==========================
while true
    % Determine step bound
    if isVectorTspan
        if nextOut > numel(T), break; end
        h = min([h, MaxStep, T(nextOut) - t]);
    else
        h = min([h, MaxStep, tspan(2) - t]);
    end
    if h <= eps(max(1,abs(t)))
        warning('ode45_cell:StepUnderflow','Step size underflow, terminating.');
        break;
    end

    % ---------- Stage evaluations k1..k6 ------------------------------
    K1 = odefun(t,                y0);
    K2 = odefun(t+C(2)*h, add (y0, h*A(2,1), K1));
    K3 = odefun(t+C(3)*h, add3(y0, h*A(3,1), K1, h*A(3,2), K2));
    K4 = odefun(t+C(4)*h, add4(y0, h*A(4,1), K1, h*A(4,2), K2, h*A(4,3), K3));
    K5 = odefun(t+C(5)*h, add5(y0, h*A(5,1), K1, h*A(5,2), K2, h*A(5,3), K3, h*A(5,4), K4));
    K6 = odefun(t+C(6)*h, add6(y0, h*A(6,1), K1, h*A(6,2), K2, h*A(6,3), K3, h*A(6,4), K4, h*A(6,5), K5));

    % Fifth‑ & fourth‑order solutions
    y5 = add7(y0, h*B5(1),K1, h*B5(3),K3, h*B5(4),K4, h*B5(5),K5, h*B5(6),K6);
    y4 = add7(y0, h*B4(1),K1, h*B4(3),K3, h*B4(4),K4, h*B4(5),K5, h*B4(6),K6, h*B4(7),K6);

    % Error estimate (real, positive)
    err = errorNorm(y4,y5,AbsTol,RelTol);

    if err <= 1          % ---------- Accept step --------------------
        t  = t + h;      % advance time
        y0 = y5;         % advance state

        if isVectorTspan
            if abs(t - T(nextOut)) <= 1e-12*max(1,abs(t))
                for i = 1:numel(Y), Y{i}(nextOut,:) = y0{i}(:).'; end
                if ~isempty(OutputFcn)
                    if feval(OutputFcn,t,y0,''), break; end
                end
                nextOut = nextOut + 1;
            end
            if t >= T(end) - eps(max(1,abs(T(end)))), break; end
        else
            % Store every accepted step
            storeIdx = storeIdx + 1;
            if storeIdx > numel(T)
                % grow storage
                T  = [T; zeros(chunk,1)];
                Y  = cellfun(@(Yi) cat(1, Yi, zeros([chunk, size(Yi,2:size(Yi,2))])), Y,'UniformOutput',false);
            end
            T(storeIdx) = t;
            for i = 1:numel(Y), Y{i}(storeIdx,:) = y0{i}(:).'; end
            if ~isempty(OutputFcn)
                if feval(OutputFcn,t,y0,''), break; end
            end
            if t >= tspan(2) - eps(max(1,abs(tspan(2)))), break; end
        end
    end

    % ---------- Update step size --------------------------------------
    h = 0.9 * h * max(0.2, min(5, err^(-1/(METHOD_ORDER+1))));
end

% Final OutputFcn call
if ~isempty(OutputFcn)
    feval(OutputFcn,t,y0,'done');
end

% Trim pre‑allocation if needed
if ~isVectorTspan
    T = T(1:storeIdx);
    Y = cellfun(@(Yi) Yi(1:storeIdx,:), Y,'UniformOutput',false);
end

end %──────────────────────────────────────────────────────────────────────

%% ======================= Local helper functions ========================
function val = getOpt(s,field,default)
    if isfield(s,field), val = s.(field); else, val = default; end
end

% ----- vectorised cell arithmetic helpers ------------------------------
function y = add(y0,a1,k1)
    y = cellfun(@(y_,k_) y_ + a1*k_, y0,k1,'UniformOutput',false);
end
function y = add3(y0,a1,k1,a2,k2)
    y = cellfun(@(y_,k1_,k2_) y_ + a1*k1_ + a2*k2_, y0,k1,k2,'UniformOutput',false);
end
function y = add4(y0,a1,k1,a2,k2,a3,k3)
    y = cellfun(@(y_,k1_,k2_,k3_) y_ + a1*k1_ + a2*k2_ + a3*k3_, y0,k1,k2,k3,'UniformOutput',false);
end
function y = add5(y0,a1,k1,a2,k2,a3,k3,a4,k4)
    y = cellfun(@(y_,k1_,k2_,k3_,k4_) y_ + a1*k1_ + a2*k2_ + a3*k3_ + a4*k4_, y0,k1,k2,k3,k4,'UniformOutput',false);
end
function y = add6(y0,a1,k1,a2,k2,a3,k3,a4,k4,a5,k5)
    y = cellfun(@(y_,k1_,k2_,k3_,k4_,k5_) y_ + a1*k1_ + a2*k2_ + a3*k3_ + a4*k4_ + a5*k5_, y0,k1,k2,k3,k4,k5,'UniformOutput',false);
end
function y = add7(y0,a1,k1,a3,k3,a4,k4,a5,k5,a6,k6,a7,k7)
    % last pair optional
    if nargin < 15, a7 = 0; k7 = k6; end
    y = cellfun(@(y_,k1_,k3_,k4_,k5_,k6_,k7_) y_ + a1*k1_ + a3*k3_ + a4*k4_ + a5*k5_ + a6*k6_ + a7*k7_,y0,k1,k3,k4,k5,k6,k7,'UniformOutput',false);
end

% ----- error norm: real RMS of magnitudes ------------------------------
function err = errorNorm(y4,y5,AbsTol,RelTol)
    scale = cellfun(@(a,b,at) at + max(abs(a),abs(b))*RelTol, y4,y5,AbsTol,'UniformOutput',false);
    E     = cellfun(@(a,b,s) abs(b - a)./s, y4,y5,scale,'UniformOutput',false);
    err   = sqrt(sum(cellfun(@(e) sum(e(:).^2), E)) / numel(E{1}(:)));
end

% ----- choose initial step --------------------------------------------
function h = initialStep(odefun,t0,y0,RelTol,AbsTol,hUser,methodOrder)
    if ~isempty(hUser), h = hUser; return; end
    f0 = odefun(t0,y0);
    scale = cellfun(@(y,at) at + abs(y)*RelTol, y0,AbsTol,'UniformOutput',false);
    d0 = sqrt(sum(cellfun(@(f,s) sum(abs(f(:)).^2 ./ s(:).^2), f0,scale)) / numel(scale{1}(:)));
    if d0 == 0, h = 1e-3; else, h = 0.01 * d0^(-1/(methodOrder+1)); end
    h = max(h,1e-6);
end
