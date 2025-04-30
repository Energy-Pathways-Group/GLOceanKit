function [T,Y] = ode45_cell(odefun,tspan,y0,options)
%ODE45_CELL  Dormand–Prince (4,5) integrator that works with CELL arrays.
%
%   This is a drop‑in replacement for MATLAB® ODE45 *except* that the state
%   Y, the initial condition Y0, and absolute tolerances may all be *cell
%   arrays*.  Each cell element can contain a scalar, vector, matrix or
%   N‑D numeric array.
%
%   [T,Y] = ODE45_CELL(ODEFUN,TSPAN,Y0,OPTIONS)
%       ODEFUN : @(t,y) -> dy/dt  (returns cell same size as Y0)
%       TSPAN  : Either [t0 tf] or a vector [t0 t1 … tf]
%       Y0     : Cell array initial state
%       OPTIONS: struct with fields (all optional)
%                   RelTol         scalar   (default 1e‑3)
%                   AbsTol | absTolerance   scalar or cell (default 1e‑6)
%                   InitialStep    scalar   first step size
%                   MaxStep        scalar   upper bound on step
%                   OutputFcn      @(t,y,state) user callback
%                                   state = 'init','', or 'done'
%
%   If TSPAN is a vector with >2 elements, the solver returns the solution
%   exactly at those times *and* calls OutputFcn at each of those points.
%
%   Example (scalar + vector stored in a cell)
%   -----------------------------------------
%   f = @(t,y){ sin(t)*y{1}; [0 1;-1 0]*y{2} };
%   [t,y] = ode45_cell(f,0:0.2:10,{1;[1;0]});
%   plot(t,squeeze(y{1})), title('Cell element 1')
%
%   ---------------------------------------------------------------------

% Dormand–Prince coefficients -------------------------------------------
A = [0 0 0 0 0 0;
     1/5 0 0 0 0 0;
     3/40 9/40 0 0 0 0;
     44/45 -56/15 32/9 0 0 0;
     19372/6561 -25360/2187 64448/6561 -212/729 0 0;
     9017/3168 -355/33 46732/5247 49/176 -5103/18656 0];
B5 = [35/384 0 500/1113 125/192 -2187/6784 11/84 0];           % 5th
B4 = [5179/57600 0 7571/16695 393/640 -92097/339200 187/2100 1/40]; % 4th
C  = [0;1/5;3/10;4/5;8/9;1];
order = 5;                        % method order for step‑control

%% Input parsing ---------------------------------------------------------
if nargin<4, options = struct(); end
if ~iscell(y0), error('y0 must be a cell array'); end
if ~isa(odefun,'function_handle'), error('odefun must be a function handle'); end

isTspanVector = numel(tspan) > 2;
if ~isTspanVector && numel(tspan)~=2
    error('tspan must be [t0 tf] or a monotone vector');
end

RelTol      = getOpt(options,'RelTol',1e-3);
AbsTolOpt   = getOpt(options,'AbsTol',[]);
AbsTolAlt   = getOpt(options,'absTolerance',[]);
if isempty(AbsTolOpt), AbsTolOpt = AbsTolAlt; end
if isempty(AbsTolOpt), AbsTolOpt = 1e-6; end
InitialStep = getOpt(options,'InitialStep',[]);
MaxStep     = getOpt(options,'MaxStep',max(tspan)-min(tspan));
OutputFcn   = getOpt(options,'OutputFcn',[]);

% Expand AbsTol to cell
if ~iscell(AbsTolOpt)
    AbsTol = cellfun(@(y) repmat(AbsTolOpt,size(y)), y0,'UniformOutput',false);
else
    if ~isequal(size(AbsTolOpt),size(y0))
        error('AbsTol/absTolerance cell array must match y0');
    end
    AbsTol = AbsTolOpt;
end

%% Output storage --------------------------------------------------------
if isTspanVector                     % user‑requested output times
    T = tspan(:);
    nOut = numel(T);
    Y = cellfun(@(y) zeros([nOut,size(y)]), y0,'UniformOutput',false);
    % fill first value
    for k = 1:numel(Y), Y{k}(1,:) = y0{k}(:).'; end
    nextOut = 2;
else                                  % store every accepted step
    alloc = 1000;
    T = zeros(alloc,1);
    Y = cellfun(@(y) zeros([alloc,size(y)]), y0,'UniformOutput',false);
    T(1) = tspan(1); for k=1:numel(Y), Y{k}(1,:) = y0{k}(:).'; end
    idxIn = 1;                         % last stored index
end

%% Initial step size -----------------------------------------------------
h = initialStep(odefun,tspan(1),y0,RelTol,AbsTol,InitialStep);

%% OutputFcn: init -------------------------------------------------------
if ~isempty(OutputFcn)
    if feval(OutputFcn,tspan(1),y0,'init'), return; end
end

t = tspan(1);

%% Main loop -------------------------------------------------------------
while true
    if isTspanVector
        if nextOut > numel(T), break; end
        tEnd = T(end);
        h = min([h, MaxStep, T(nextOut)-t]);
    else
        tEnd = tspan(2);
        h = min([h, MaxStep, tEnd - t]);
    end
    if h <= eps(max(1,abs(t)))
        warning('Step size underflow. Integration terminated early.');
        break;
    end

    %-------------------------------- Stage evaluations k1..k6 ----------
    K1 = feval(odefun,t,y0);
    K2 = feval(odefun,t+C(2)*h, add(y0,h*A(2,1),K1));
    K3 = feval(odefun,t+C(3)*h, add3(y0,h*A(3,1),K1,h*A(3,2),K2));
    K4 = feval(odefun,t+C(4)*h, add4(y0,h*A(4,1),K1,h*A(4,2),K2,h*A(4,3),K3));
    K5 = feval(odefun,t+C(5)*h, add5(y0,h*A(5,1),K1,h*A(5,2),K2,h*A(5,3),K3,h*A(5,4),K4));
    K6 = feval(odefun,t+C(6)*h, add6(y0,h*A(6,1),K1,h*A(6,2),K2,h*A(6,3),K3,h*A(6,4),K4,h*A(6,5),K5));

    y5 = add7(y0,h*B5(1),K1,h*B5(3),K3,h*B5(4),K4,h*B5(5),K5,h*B5(6),K6);
    y4 = add7(y0,h*B4(1),K1,h*B4(3),K3,h*B4(4),K4,h*B4(5),K5,h*B4(6),K6,h*B4(7),K6);

    err = errorNorm(y4,y5,AbsTol,RelTol);

    if err <= 1          % Accept step -------------------------------
        t  = t + h;
        y0 = y5;

        if isTspanVector
            if abs(t - T(nextOut)) <= 1e-12*max(1,abs(t))
                % hit requested output point exactly
                for k = 1:numel(Y), Y{k}(nextOut,:) = y0{k}(:).'; end
                if ~isempty(OutputFcn)
                    if feval(OutputFcn,t,y0,''), break; end
                end
                nextOut = nextOut + 1;
            end
            if t >= T(end) - eps(max(1,abs(T(end))))
                break;   % reached final user time
            end
        else
            idxIn = idxIn + 1;
            if idxIn > numel(T)
                T = [T; zeros(1000,1)];
                Y = cellfun(@(Yi) cat(1,Yi,zeros([1000,size(Yi,2:size(Yi,2))])), Y,'UniformOutput',false);
            end
            T(idxIn) = t;
            for k = 1:numel(Y), Y{k}(idxIn,:) = y0{k}(:).'; end
            if ~isempty(OutputFcn)
                if feval(OutputFcn,t,y0,''), break; end
            end
            if t >= tEnd - eps(max(1,abs(tEnd)))
                break;
            end
        end
    end

    % Step‑size update ---------------------------------------------------
    h = 0.9*h*max(0.2,min(5, err^(-1/(order+1))));
end

%% OutputFcn: done -------------------------------------------------------
if ~isempty(OutputFcn)
    feval(OutputFcn,t,y0,'done');
end

%% Trim dynamic storage --------------------------------------------------
if isTspanVector
    % already sized correctly
else
    T = T(1:idxIn);
    Y = cellfun(@(Yi) Yi(1:idxIn,:), Y,'UniformOutput',false);
end

end % ode45_cell

%% ────────────────── Helper functions ───────────────────────────────────
function val = getOpt(s,f,d), if isfield(s,f), val = s.(f); else, val = d; end, end

% Vectorised cell additions ------------------------------------------------
function y = add(y0,a1,k1)
    y = cellfun(@(y,k) y + a1*k, y0,k1,'UniformOutput',false);
end
function y = add3(y0,a1,k1,a2,k2)
    y = cellfun(@(y,k1i,k2i) y + a1*k1i + a2*k2i, y0,k1,k2,'UniformOutput',false);
end
function y = add4(y0,a1,k1,a2,k2,a3,k3)
    y = cellfun(@(y,k1i,k2i,k3i) y + a1*k1i + a2*k2i + a3*k3i, y0,k1,k2,k3,'UniformOutput',false);
end
function y = add5(y0,a1,k1,a2,k2,a3,k3,a4,k4)
    y = cellfun(@(y,k1i,k2i,k3i,k4i) y + a1*k1i + a2*k2i + a3*k3i + a4*k4i, y0,k1,k2,k3,k4,'UniformOutput',false);
end
function y = add6(y0,a1,k1,a2,k2,a3,k3,a4,k4,a5,k5)
    y = cellfun(@(y,k1i,k2i,k3i,k4i,k5i) y + a1*k1i + a2*k2i + a3*k3i + a4*k4i + a5*k5i, y0,k1,k2,k3,k4,k5,'UniformOutput',false);
end
function y = add7(y0,a1,k1,a3,k3,a4,k4,a5,k5,a6,k6,a7,k7)
    if nargin<15, a7=0; k7=k6; end
    y = cellfun(@(y,k1i,k3i,k4i,k5i,k6i,k7i) y + a1*k1i + a3*k3i + a4*k4i + a5*k5i + a6*k6i + a7*k7i, y0,k1,k3,k4,k5,k6,k7,'UniformOutput',false);
end

% Error norm --------------------------------------------------------------
function err = errorNorm(y4,y5,AbsTol,RelTol)
    denom = cellfun(@(y,absi) absi + max(abs(y),abs(y5))*RelTol, y5,AbsTol,'UniformOutput',false);
    E = cellfun(@(s4,s5,d) (s5-s4)./d, y4,y5,denom,'UniformOutput',false);
    err = sqrt(sum(cellfun(@(e) sum(e(:).^2),E))/numel(E{1}(:)));
end

% Initial step heuristic --------------------------------------------------
function h = initialStep(odefun,t0,y0,RelTol,AbsTol,hUser)
    if ~isempty(hUser), h = hUser; return; end
    f0 = odefun(t0,y0);
    scale = cellfun(@(y,a) a + abs(y)*RelTol, y0,AbsTol,'UniformOutput',false);
    d0 = sqrt(sum(cellfun(@(f,s) sum((f(:)./s(:)).^2), f0,scale))/numel(scale{1}(:)));
    h = 0.01 * max(1e-6, 1/d0^(1/(order+1)));
end
