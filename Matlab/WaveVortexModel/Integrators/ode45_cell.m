function ode45_cell(odefun,tspan,y0,options)
%ODE45_CELL  Dormand–Prince (4,5) integrator for CELL-array states.
%
%   ode45_cell(ODEFUN,TSPAN,Y0) integrates y' = f(t,y) where Y is a cell
%   array.  The function **does not return** any solution vectors; you must
%   use an OutputFcn to access the evolving state.  Complex numbers are
%   supported.
%
%   INPUTS
%   ------
%   ODEFUN  @(t,y) -> dy/dt (cell array matching Y0)
%   TSPAN   [t0 tf] or a vector of output times
%   Y0      cell array initial condition
%   OPTIONS struct (all fields optional)
%              RelTol        (1e-3)
%              AbsTol        (1e-6 scalar or cell)
%              InitialStep   ([] automatic)
%              MaxStep       (tf-t0)
%              OutputFcn(t,y,state)  state='init','', or 'done'
% ----------------------------------------------------------------------

%% Butcher tableau (Dormand–Prince 4,5) ---------------------------------
A = [0 0 0 0 0 0;
     1/5 0 0 0 0 0;
     3/40 9/40 0 0 0 0;
     44/45 -56/15 32/9 0 0 0;
     19372/6561 -25360/2187 64448/6561 -212/729 0 0;
     9017/3168 -355/33 46732/5247 49/176 -5103/18656 0];
B5 = [35/384 0 500/1113 125/192 -2187/6784 11/84 0];
B4 = [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40];
C  = [0; 1/5; 3/10; 4/5; 8/9; 1];
ORDER = 5;

%% ---------------- Parse inputs ----------------------------------------
if nargin < 4, options = struct(); end
if ~iscell(y0),      error('Y0 must be a cell array.'); end
if ~isa(odefun,'function_handle'), error('ODEFUN must be a function handle.'); end

vecTspan = numel(tspan) > 2;
if ~vecTspan && numel(tspan)~=2
    error('TSPAN must be a 2‑element interval or a monotone vector.');
end

RelTol      = getOpt(options,'RelTol',1e-3);
AbsTolIn    = getOpt(options,'AbsTol',1e-6);
InitialStep = getOpt(options,'InitialStep',[]);
MaxStep     = getOpt(options,'MaxStep',max(tspan)-min(tspan));
OutputFcn   = getOpt(options,'OutputFcn',[]);

% Expand AbsTol to cell
if ~iscell(AbsTolIn)
    AbsTol = cellfun(@(y) repmat(AbsTolIn,size(y)), y0,'UniformOutput',false);
else
    if ~isequal(size(AbsTolIn),size(y0))
        error('AbsTol cell array must match Y0.');
    end
    AbsTol = AbsTolIn;
end

%% ---------------- Initial step ----------------------------------------
h = initialStep(odefun,tspan(1),y0,RelTol,AbsTol,InitialStep,ORDER);

t = tspan(1);
if vecTspan, nextIdx = 2; end

% OutputFcn init
if ~isempty(OutputFcn)
    if feval(OutputFcn,t,y0,'init'), return; end
end

%% ---------------- Main loop -------------------------------------------
while true
    % Cap step by remaining interval / next requested time
    if vecTspan
        if nextIdx > numel(tspan), break; end
        h = min([h, MaxStep, tspan(nextIdx)-t]);
    else
        h = min([h, MaxStep, tspan(2)-t]);
    end
    if h <= eps(max(1,abs(t)))
        warning('ode45_cell:StepTooSmall','Step size underflow.');
        break;
    end

    % ----- Stages -------------------------------------------------------
    k1 = odefun(t,                y0);
    k2 = odefun(t+C(2)*h, add (y0,h*A(2,1),k1));
    k3 = odefun(t+C(3)*h, add3(y0,h*A(3,1),k1,h*A(3,2),k2));
    k4 = odefun(t+C(4)*h, add4(y0,h*A(4,1),k1,h*A(4,2),k2,h*A(4,3),k3));
    k5 = odefun(t+C(5)*h, add5(y0,h*A(5,1),k1,h*A(5,2),k2,h*A(5,3),k3,h*A(5,4),k4));
    k6 = odefun(t+C(6)*h, add6(y0,h*A(6,1),k1,h*A(6,2),k2,h*A(6,3),k3,h*A(6,4),k4,h*A(6,5),k5));

    y5 = add7(y0,h*B5(1),k1,h*B5(3),k3,h*B5(4),k4,h*B5(5),k5,h*B5(6),k6);
    y4 = add7(y0,h*B4(1),k1,h*B4(3),k3,h*B4(4),k4,h*B4(5),k5,h*B4(6),k6,h*B4(7),k6);

    err = errorNorm(y4,y5,AbsTol,RelTol);
    if err <= 1  % Accept step
        t  = t + h;
        y0 = y5;

        atOutput = (~vecTspan) || (abs(t - tspan(nextIdx)) <= 1e-12*max(1,abs(t)));
        if atOutput && ~isempty(OutputFcn)
            if feval(OutputFcn,t,y0,''), break; end
        end
        if vecTspan && atOutput
            nextIdx = nextIdx + 1;
            if t >= tspan(end) - eps(max(1,abs(tspan(end)))), break; end
        elseif ~vecTspan && t >= tspan(2) - eps(max(1,abs(tspan(2)))), break; end
    end

    % Update step size
    % h = 0.9*h*max(0.2,min(5, err^(-1/(ORDER+1))));
    safety   = 0.9;        % like MATLAB's ode45
    maxGrow  = 5;          % h ≤ 5·h_old  (MATLAB uses ~5)
    minShrink= 0.1;        % h ≥ 0.1·h_old (MATLAB uses 0.1)

    h = safety * h * max(minShrink, min(maxGrow, err^(-1/5)));
end

% OutputFcn done
if ~isempty(OutputFcn)
    feval(OutputFcn,t,y0,'done');
end

end % =================== primary function end ===========================

%% ----------------- Helper functions -----------------------------------
function val = getOpt(s,f,d)
    if isfield(s,f), val = s.(f); else, val = d; end
end

% Vectorised cell arithmetic
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
    if nargin<15, a7 = 0; k7 = k6; end
    y = cellfun(@(y_,k1_,k3_,k4_,k5_,k6_,k7_) y_ + a1*k1_ + a3*k3_ + a4*k4_ + a5*k5_ + a6*k6_ + a7*k7_,...
        y0,k1,k3,k4,k5,k6,k7,'UniformOutput',false);
end

% Error norm (max of magnitudes)
function err = errorNorm(y4,y5,AbsTol,RelTol)
    % scale = cellfun(@(a,at)  max(abs(a)*RelTol,at), y5,AbsTol,'UniformOutput',false);
    scale = cellfun(@(a,b,at) at + max(abs(a),abs(b))*RelTol, y4, y5, AbsTol, 'UniformOutput', false);
    E     = cellfun(@(a,b,s) abs(b-a)./s, y4,y5,scale,'UniformOutput',false);
    err   = max(cellfun(@(e) max(e(:)), E));
end

% Initial step size heuristic
function h = initialStep(odefun,t0,y0,RelTol,AbsTol,hUser,order)
    if ~isempty(hUser), h = hUser; return; end
    f0 = odefun(t0,y0);
    scale = cellfun(@(y,a) a + abs(y)*RelTol, y0,AbsTol,'UniformOutput',false);
    d0 = sqrt(sum(cellfun(@(f,s) sum(abs(f(:)).^2 ./ s(:).^2), f0,scale)) / numel(scale{1}(:)));
    if d0 == 0, h = 1e-3; else, h = 0.01 * d0^(-1/(order+1)); end
    h = max(h,1e-6);
end
