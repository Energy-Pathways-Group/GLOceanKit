function ode45_cell(odefun, tspan, y0, options)
%
% 2025-05-01
% I tried using ChatGPT to create this and it was a mixed bag. I believe it
% incorrectly did not compute k7, and even upon my insistance, it said no,
% I was wrong. Adding the k7 stage and making it function as I read the
% Butcher table, everything now matches the time step in Matlab (previously
% it was about 1/10 the time step). Thus I think the Dormand-Prince is
% basically working. the Bogacki-Shampine, however, is not. I will come
% back to this later.
%
%ODE45_CELL   Adaptive RK(4,5) solver for *cell-array* states (no output vectors).
%   • OPTIONS.Method      = 'DP' (Dormand‑Prince, default)  | 'BS' (Bogacki‑Shampine)
%   • OPTIONS.ErrorNorm   = 'rms' (MATLAB default) | 'max'
%   • OPTIONS.RelTol      = 1e-3   (scalar)
%   • OPTIONS.AbsTol      = 1e-6   (scalar or cell the size of Y0)
%   • OPTIONS.InitialStep = []     (solver chooses)
%   • OPTIONS.MaxStep     = tf–t0  (scalar)
%   • OPTIONS.OutputFcn(t,y,flag)  where flag = 'init' | '' | 'done'
%
%   The solver stores no solution: use an OutputFcn to capture results.
% -----------------------------------------------------------------------------

%% ----------------------- Option handling ----------------------------------
if nargin < 4
    options = struct();
end
Method     = getOpt(options, 'Method',     'DP');     % 'DP' | 'BS'
ErrChoice  = getOpt(options, 'ErrorNorm',  'rms');    % 'rms'| 'max'
RelTol     = getOpt(options, 'RelTol',      1e-3);
AbsTol_in  = getOpt(options, 'AbsTol',      1e-6);
InitialStep= getOpt(options, 'InitialStep', []   );
MaxStep    = getOpt(options, 'MaxStep',     max(tspan)-min(tspan));
OutputFcn  = getOpt(options, 'OutputFcn',   []   );

%% ----------------------- Checks -------------------------------------------
assert(iscell(y0), 'Y0 must be a cell array');
assert(isa(odefun,'function_handle'),'ODEFUN must be a function handle');
vecGrid = numel(tspan) > 2;
if ~vecGrid
    assert(numel(tspan)==2,'TSPAN must be [t0 tf] or a vector');
end
if ~iscell(AbsTol_in)
    AbsTol = cellfun(@(y) repmat(AbsTol_in, size(y)), y0, 'uni', 0);
else
    assert(isequal(size(AbsTol_in), size(y0)), 'AbsTol cell size mismatch');
    AbsTol = AbsTol_in;
end

%% ----------------------- Coefficient tables --------------------------------
if strcmpi(Method,'DP')
    C = [ 0;
          1/5;
          3/10;
          4/5;
          8/9;
          1 ];

    A = [ 0           0           0           0           0           0;
          1/5         0           0           0           0           0;
          3/40        9/40        0           0           0           0;
          44/45     -56/15       32/9         0           0           0;
       19372/6561 -25360/2187 64448/6561   -212/729       0           0;
        9017/3168    -355/33   46732/5247     49/176  -5103/18656     0 ];

    B5 = [ 35/384   0   500/1113   125/192   -2187/6784    11/84    0 ];
    B4 = [ 5179/57600   0   7571/16695   393/640  -92097/339200  187/2100  1/40 ];
else                               % Bogacki–Shampine 5(4)
    C = [ 0;
          2/9;
          1/3;
          3/4;
          1;
          5/6 ];

    A = [ 0           0           0           0         0        0;
          2/9         0           0           0         0        0;
          1/12        1/4         0           0         0        0;
          69/128    -243/128     135/64       0         0        0;
         -17/12       27/4      -27/5       16/15       0        0;
          65/432     -5/16       13/16       4/27     5/144      0 ];

    B5 = [ 47/450   0   12/25   32/225   1/30   6/25   0 ];
    B4 = [  6/25   0  -40/27   82/9    32/54   4/15   0 ];
end

%% ----------------------- Initial step -------------------------------------
if isempty(InitialStep)
    h = initStep(odefun, tspan(1), y0, RelTol, AbsTol);
else
    h = InitialStep;
end
h = max(h, 1e-6);

%% ----------------------- Control parameters --------------------------------
ORDER = 5;           % order of high solution for exponent
safety = 0.8;        % MATLAB default
minShrink = 0.2; maxGrow = 10;

useRMS = strcmpi(ErrChoice,'rms');

prev_t  = tspan(1);
prev_y  = y0;
prev_k  = odefun(prev_t, y0);

t = prev_t;
if vecGrid, nextIdx = 2; end
if ~isempty(OutputFcn)
    if OutputFcn(t, y0, 'init'); return; end
end

%% ======================= MAIN LOOP =========================================
while true
    h = min(h, MaxStep);
    if ~vecGrid
        h = min(h, tspan(2) - t);
    end
    if h <= eps(max(1, abs(t)))
        warning('ode45_cell:StepTooSmall', 'Step size too small');
        break
    end

    % ---- Stage evaluations -------------------------------------------
    k1 = prev_k;
    k2 = odefun(t + C(2)*h, add   (prev_y, h*A(2,1), k1));
    k3 = odefun(t + C(3)*h, add3  (prev_y, h*A(3,1), k1, h*A(3,2), k2));
    k4 = odefun(t + C(4)*h, add4  (prev_y, h*A(4,1), k1, h*A(4,2), k2, h*A(4,3), k3));
    k5 = odefun(t + C(5)*h, add5  (prev_y, h*A(5,1), k1, h*A(5,2), k2, h*A(5,3), k3, h*A(5,4), k4));
    k6 = odefun(t + C(6)*h, add6  (prev_y, h*A(6,1), k1, h*A(6,2), k2, h*A(6,3), k3, h*A(6,4), k4, h*A(6,5), k5));

    y5 = add7(prev_y, h*B5(1), k1, h*B5(3), k3, h*B5(4), k4, h*B5(5), k5, h*B5(6), k6);
    k7 = odefun(t + C(6)*h, y5);
    y4 = add7(prev_y, h*B4(1), k1, h*B4(3), k3, h*B4(4), k4, h*B4(5), k5, h*B4(6), k6, h*B4(7), k7);

    if useRMS
        err = rmsError(y4, y5, AbsTol, RelTol);
    else
        err = maxError(y4, y5, AbsTol, RelTol);
    end

    if err <= 1
        % -------- accept step ----------------------------------------
        new_t = t + h;
        new_y = y5;
        new_k = k7;

        if vecGrid && ~isempty(OutputFcn)
            while nextIdx <= numel(tspan) && tspan(nextIdx) <= new_t + eps(max(1, abs(new_t)))
                tau   = (tspan(nextIdx) - t) / h;
                y_out = hermiteInterp(prev_y, new_y, prev_k, new_k, tau, h);
                if OutputFcn(tspan(nextIdx), y_out, '')
                    return
                end
                nextIdx = nextIdx + 1;
            end
        elseif ~vecGrid && ~isempty(OutputFcn)
            if OutputFcn(new_t, new_y, '')
                break
            end
        end

        if (~vecGrid && new_t >= tspan(2) - eps(max(1, abs(tspan(2))))) || ...
           ( vecGrid && nextIdx > numel(tspan) )
            break
        end

        prev_t = new_t;
        prev_y = new_y;

        % Select correct derivative at the end of the step
        if strcmpi(Method,'DP')
            % Dormand–Prince: stage k6 is evaluated at c = 1
            prev_k = k7;
        else
            % Bogacki–Shampine: stage k5 is evaluated at c = 1
            prev_k = k5;
        end
        t = new_t;
    % else
    %     fprintf('error too large\n');
    end

    % ---- Step-size controller --------------------------------------
    h = safety * h * err^(-1/ORDER);
    h = max(minShrink * h, min(maxGrow * h, h));
end

if ~isempty(OutputFcn)
    OutputFcn(t, prev_y, 'done');
end
end %==================== MAIN FUNCTION END ================================

%% ======================================================================
%% Dense-output Hermite interpolant ----------------------------------------
function y = hermiteInterp(y0, y1, k0, k1, theta, h)
    % 4th‑order cubic Hermite interpolation between two accepted steps
    h00 =  2*theta^3 - 3*theta^2 + 1;
    h10 =      theta^3 - 2*theta^2 + theta;
    h01 = -2*theta^3 + 3*theta^2;
    h11 =      theta^3 -     theta^2;

    y = cellfun(@(y0i, y1i, k0i, k1i) ...
                h00 .* y0i + h10 .* h .* k0i + h01 .* y1i + h11 .* h .* k1i, ...
                y0,  y1,  k0,  k1, 'UniformOutput', false);
end

%% ----------------------------------------------------------------------
%% Error norms
function err = rmsError(y4, y5, AbsTol, RelTol)
    % Weighted root‑mean‑square local error (MATLAB style)
    scale = cellfun(@(a, b, at) at + max(abs(a), abs(b)) .* RelTol, y4, y5, AbsTol, 'UniformOutput', false);
    sq = cellfun(@(a, b, s) (abs(b - a) ./ s).^2, y4, y5, scale, 'UniformOutput', false);
    num = sum(cellfun(@(m) sum(m(:)), sq));
    den = sum(cellfun(@numel, sq));
    err = sqrt(num / den);
end

function err = maxError(y4, y5, AbsTol, RelTol)
    % Infinity‑norm local error (conservative)
    scale = cellfun(@(a, b, at) at + max(abs(a), abs(b)) .* RelTol, ...
                    y4, y5, AbsTol, 'UniformOutput', false);
    e = cellfun(@(a, b, s) abs(b - a) ./ s, y4, y5, scale, 'UniformOutput', false);
    err = max(cellfun(@(ei) max(ei(:)), e));
end

%% ----------------------------------------------------------------------
%% Initial step heuristic
function h = initStep(odefun, t0, y0, RelTol, AbsTol)
    f0 = odefun(t0, y0);
    scale = cellfun(@(y, a) a + abs(y) .* RelTol, y0, AbsTol, 'UniformOutput', false);
    d0 = max(cellfun(@(f, s) max(abs(f(:) ./ s(:))), f0, scale));

    if d0 == 0
        h = 1e-3;
    else
        h = 0.16 / d0;   % Shampine recommendation
    end
end

%% ----------------------------------------------------------------------
%% Option helper (accepts empty field as missing)
function v = getOpt(s, field, default)
    if isfield(s, field)
        if isempty(s.(field))
            v = default;
        else
            v = s.(field);
        end
    else
        v = default;
    end
end

%% ----------------------------------------------------------------------
%% Cell‑wise arithmetic helpers (element‑wise scalar coeffs)
function y = add(y0, a1, k1)
    y = cellfun(@(y_, k_) y_ + k_ .* a1, y0, k1, 'UniformOutput', false);
end
function y = add3(y0, a1, k1, a2, k2)
    y = cellfun(@(y_, k1_, k2_) y_ + k1_ .* a1 + k2_ .* a2, ...
                y0, k1, k2, 'UniformOutput', false);
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

%% ========================= End of file =================================
