function ode45_cell(odefun, tspan, y0, options)
%ODE45_CELL  Dormand–Prince (4,5) integrator for CELL‑array states
%            (no solution arrays are returned).  Use an OutputFcn to
%            receive results.  Supports complex data and provides a
%            4‑th‑order Hermite dense‑output interpolant when TSPAN is a
%            vector.
%
%   ode45_cell(ODEFUN, TSPAN, Y0)
%   ode45_cell(ODEFUN, TSPAN, Y0, OPTIONS)
%
%   Required
%   --------
%   ODEFUN   @(t, y) cell ‑ derivative function (same cell layout as Y0)
%   TSPAN    [t0 tf] or a monotone vector of output times
%   Y0       cell array initial condition
%
%   OPTIONS fields (all optional)
%   --------------------------------
%   RelTol        — scalar (1e‑3)
%   AbsTol        — scalar or matching cell (1e‑6)
%   InitialStep   — scalar (automatic if empty)
%   MaxStep       — scalar (tf‑t0)
%   OutputFcn(t, y, state) where state = 'init' | '' | 'done'
% -------------------------------------------------------------------------

%% Butcher tableau: Dormand–Prince 5(4) ----------------------------------
A = [  0          0          0          0          0          0 ;
       1/5        0          0          0          0          0 ;
       3/40       9/40       0          0          0          0 ;
       44/45    -56/15      32/9        0          0          0 ;
    19372/6561 -25360/2187 64448/6561  -212/729     0          0 ;
     9017/3168   -355/33   46732/5247   49/176   -5103/18656   0 ];

B5 = [ 35/384      0    500/1113   125/192   -2187/6784   11/84     0 ];
B4 = [ 5179/57600  0    7571/16695  393/640   -92097/339200 187/2100 1/40 ];
C  = [ 0 ; 1/5 ; 3/10 ; 4/5 ; 8/9 ; 1 ];

%% ---------------- Parse inputs -----------------------------------------
if nargin < 4
    options = struct();
end

assert(iscell(y0), 'Y0 must be a cell array.');
assert(isa(odefun, 'function_handle'), 'ODEFUN must be a function handle.');

vecGrid = numel(tspan) > 2;
if ~vecGrid
    assert(numel(tspan) == 2, 'TSPAN must be [t0 tf] or a vector.');
end

RelTol      = getOpt(options, 'RelTol',      1e-3 );
AbsTolInput = getOpt(options, 'AbsTol',      1e-6 );
InitialStep = getOpt(options, 'InitialStep', []   );
MaxStep     = getOpt(options, 'MaxStep',     max(tspan) - min(tspan));
OutputFcn   = getOpt(options, 'OutputFcn',   []   );

if ~iscell(AbsTolInput)
    AbsTol = cellfun(@(y) repmat(AbsTolInput, size(y)), y0, 'UniformOutput', false);
else
    assert(isequal(size(AbsTolInput), size(y0)), 'AbsTol cell must match Y0');
    AbsTol = AbsTolInput;
end

%% ---------------- Initialisation ---------------------------------------
if isempty(InitialStep)
    h = initStep(odefun, tspan(1), y0, RelTol, AbsTol);
else
    h = InitialStep;
end

h = max(h, 1e-6);

t = tspan(1);
if vecGrid
    nextIdx = 2;
end

if ~isempty(OutputFcn)
    stop = OutputFcn(t, y0, 'init');
    if stop
        return;
    end
end

prev_t = t;
prev_y = y0;
prev_k = odefun(t, y0);

%% ---------------- Main integration loop --------------------------------
while true
    % Limit step by MaxStep and, for two‑point span, by end time
    h = min(h, MaxStep);
    if ~vecGrid
        h = min(h, tspan(2) - t);
    end
    if h <= eps(max(1, abs(t)))
        warning('ode45_cell:StepUnderflow', 'Step size too small.');
        break;
    end

    % ----- Stage evaluations -------------------------------------------
    k1 = prev_k;
    k2 = odefun(t + C(2)*h, add   (y0, h*A(2,1), k1));
    k3 = odefun(t + C(3)*h, add3  (y0, h*A(3,1), k1, h*A(3,2), k2));
    k4 = odefun(t + C(4)*h, add4  (y0, h*A(4,1), k1, h*A(4,2), k2, h*A(4,3), k3));
    k5 = odefun(t + C(5)*h, add5  (y0, h*A(5,1), k1, h*A(5,2), k2, h*A(5,3), k3, h*A(5,4), k4));
    k6 = odefun(t + C(6)*h, add6  (y0, h*A(6,1), k1, h*A(6,2), k2, h*A(6,3), k3, h*A(6,4), k4, h*A(6,5), k5));

    y5 = add7(y0, h*B5(1), k1, h*B5(3), k3, h*B5(4), k4, h*B5(5), k5, h*B5(6), k6);
    y4 = add7(y0, h*B4(1), k1, h*B4(3), k3, h*B4(4), k4, h*B4(5), k5, h*B4(6), k6, h*B4(7), k6);

    % err = maxError(y4, y5, AbsTol, RelTol);
    err = rmsError(y4, y5, AbsTol, RelTol);

    if err <= 1
        % ---------------- Accept step ----------------------------------
        new_t = t + h;
        new_y = y5;
        new_k = k6;

        % Dense output for requested grid points
        if vecGrid && ~isempty(OutputFcn)
            while nextIdx <= numel(tspan) && tspan(nextIdx) <= new_t + eps(max(1, abs(new_t)))
                tau   = (tspan(nextIdx) - t) / h;            % fraction in [0,1]
                y_out = hermiteInterp(prev_y, new_y, prev_k, new_k, tau, h);
                stop  = OutputFcn(tspan(nextIdx), y_out, '');
                if stop
                    return;
                end
                nextIdx = nextIdx + 1;
            end
        elseif ~vecGrid && ~isempty(OutputFcn)
            stop = OutputFcn(new_t, new_y, '');
            if stop
                break;
            end
        end

        % Termination tests
        if ~vecGrid
            if new_t >= tspan(2) - eps(max(1, abs(tspan(2))))
                break;
            end
        else
            if nextIdx > numel(tspan)
                break;
            end
        end

        % Advance to next step
        prev_t = new_t;
        prev_y = new_y;
        prev_k = new_k;
        t      = new_t;
        y0     = new_y;
    end

    % Step‑size update (MATLAB‑style)
    % h = 0.8 * h * max(0.2, min(5, err^(-1/5)));

    % safety     = 0.7;   % was 0.8
    % minShrink  = 0.5;   % was 0.2   (h can shrink at most 2×)
    % maxGrow    = 2.5;   % was 10    (h can grow at most 2.5×)
    % h = safety * h * max(minShrink, min(maxGrow, err^(-1/5)));

    % --- step-size controller (inside main loop) -----------------------
    h_new = 0.8 * h * err^(-1/5);           % base formula
    h_new = max(0.5*h, min(10*h, h_new));   % conventional 0.5–10 limits
    h_new = min(1.2*h, h_new);              % <-- damped growth: at most ×1.2
    h     = h_new;
end

if ~isempty(OutputFcn)
    OutputFcn(t, y0, 'done');
end
end % ==== main function ==================================================

%% ======================================================================
%% Dense‑output: 4th‑order Hermite polynomial
function y = hermiteInterp(y0, y1, k0, k1, theta, h)
    h00 =  2*theta^3 - 3*theta^2 + 1;
    h10 =      theta^3 - 2*theta^2 + theta;
    h01 = -2*theta^3 + 3*theta^2;
    h11 =      theta^3 -     theta^2;

    y = cellfun(@(y0i, y1i, k0i, k1i) ...
                h00 .* y0i + h10 .* h .* k0i + h01 .* y1i + h11 .* h .* k1i, ...
                y0, y1, k0, k1, 'UniformOutput', false);
end

%% ======================================================================
%% Error norm (max of scaled errors)
function err = maxError(y4, y5, AbsTol, RelTol)
    scale = cellfun(@(a, b, at) at + max(abs(a), abs(b)) .* RelTol, ...
                    y4, y5, AbsTol, 'UniformOutput', false);
    e     = cellfun(@(a, b, s) abs(b - a) ./ s, y4, y5, scale, 'UniformOutput', false);
    err   = max(cellfun(@(ei) max(ei(:)), e));
end

%% ======================================================================
%% Error norm (max of scaled errors)
function err = rmsError(y4, y5, AbsTol, RelTol)
    % Weighted root-mean-square local error (matches MATLAB ode45)
    scale = cellfun(@(a, b, at) at + max(abs(a), abs(b)) .* RelTol, y4, y5, AbsTol, 'UniformOutput', false);
    sq = cellfun(@(a, b, s) (abs(b - a) ./ s).^2, y4, y5, scale, 'UniformOutput', false);

    num = sum(cellfun(@(m) sum(m(:)), sq));   % Σ error_i^2
    den = sum(cellfun(@numel, sq));           % total components N
    err = sqrt(num / den);                    % RMS error
end

%% ======================================================================
%% Initial step heuristic (simple)
function h = initStep(odefun, t0, y0, RelTol, AbsTol)
    f0 = odefun(t0, y0);
    scale = cellfun(@(y, a) a + abs(y) .* RelTol, y0, AbsTol, 'UniformOutput', false);
    d0 = max(cellfun(@(f, s) max(abs(f(:) ./ s(:))), f0, scale));

    if d0 == 0
        h = 1e-3;
    else
        h = 0.1 / d0;   % MATLAB‑like factor
    end
end

%% ======================================================================
%% Option helper
function v = getOpt(s, field, default)
    if isfield(s, field) && ~isempty(s.(field))
        v = s.(field);
    else
        v = default;
    end
end

%% ======================================================================
%% Cell‑wise arithmetic helpers (elementwise scalar coeffs)
function y = add(y0, a1, k1)
    y = cellfun(@(y_, k_) y_ + k_ .* a1, y0, k1, 'UniformOutput', false);
end
function y = add3(y0, a1, k1, a2, k2)
    y = cellfun(@(y_, k1_, k2_) y_ + k1_ .* a1 + k2_ .* a2, y0, k1, k2, 'UniformOutput', false);
end
function y = add4(y0, a1, k1, a2, k2, a3, k3)
    y = cellfun(@(y_, k1_, k2_, k3_) y_ + k1_ .* a1 + k2_ .* a2 + k3_ .* a3, ...
                y0, k1, k2, k3, 'UniformOutput', false);
end
function y = add5(y0, a1, k1, a2, k2, a3, k3, a4, k4)
    y = cellfun(@(y_, k1_, k2_, k3_, k4_) ...
                y_ + k1_ .* a1 + k2_ .* a2 + k3_ .* a3 + k4_ .* a4, ...
                y0, k1, k2, k3, k4, 'UniformOutput', false);
end
function y = add6(y0, a1, k1, a2, k2, a3, k3, a4, k4, a5, k5)
    y = cellfun(@(y_, k1_, k2_, k3_, k4_, k5_) ...
                y_ + k1_ .* a1 + k2_ .* a2 + k3_ .* a3 + k4_ .* a4 + k5_ .* a5, ...
                y0, k1, k2, k3, k4, k5, 'UniformOutput', false);
end
function y = add7(y0, a1, k1, a3, k3, a4, k4, a5, k5, a6, k6, a7, k7)
    if nargin < 13
        a7 = 0;
        k7 = k6;
    end

    y = cellfun(@(y_, k1_, k3_, k4_, k5_, k6_, k7_) ...
                y_ + k1_ .* a1 + k3_ .* a3 + k4_ .* a4 + k5_ .* a5 + k6_ .* a6 + k7_ .* a7, ...
                y0, k1, k3, k4, k5, k6, k7, 'UniformOutput', false);
end

% ======================================================================
% End of helper functions
