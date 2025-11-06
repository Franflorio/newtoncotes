function [I, out] = simpson13_selector(f, a, b, tol, opts)
% simpson13_selector  Template que elige entre apriori, halving uniforme o adaptativo.
% Inputs:
%   f     handle a funcion escalar, ej: @(x) sin(x)
%   a,b   limites con b > a
%   tol   tolerancia absoluta deseada
%   opts  (opcional) struct con campos:
%         .M4           cota de |f''''(x)| en [a,b] (>=0, finita) para modo apriori
%         .n0           n inicial para halving (default 10, se ajusta a par)
%         .kmax_auto    max refinamientos en halving global (default 20)
%         .hmin_adapt   ancho minimo en adaptativo (default 0)
%         .maxit_adapt  max subdivisiones en adaptativo (default 1e6)
%
% Outputs:
%   I     aproximacion de la integral
%   out   struct con auditoria:
%         .metodo          'apriori' | 'uniforme_auto' | 'adaptativo'
%         .h               paso final si uniforme, NaN si adaptativo
%         .n               subintervalos si uniforme, NaN si adaptativo
%         .err_est         estimacion global si uniforme, NaN si adaptativo
%         .M4_used         valor de M4 usado (o NaN)
%         .fevals_total    evaluaciones totales de f
%         .detalle         substruct con info del metodo seleccionado
%
% Notas:
%   - Apriori usa cota: |E| <= (b-a)/180 * h^4 * M4, ajusta n a par y verifica con halving 1 paso.
%   - Halving global usa err_est = |I_{2n} - I_n| / 15 (orden 4).
%   - Adaptativo usa test local 1/15 y correccion local; devuelve malla en .detalle.

    if nargin < 4
        error('Se requieren f, a, b, tol.');
    end
    if ~isa(f, 'function_handle')
        error('f debe ser function handle.');
    end
    if ~isscalar(a) || ~isscalar(b) || ~(b > a)
        error('Se requiere a < b.');
    end
    if ~isscalar(tol) || ~(tol > 0) || ~isfinite(tol)
        error('tol debe ser escalar positivo finito.');
    end

    if nargin < 5 || isempty(opts)
        opts = struct();
    end
    if ~isfield(opts, 'n0'); opts.n0 = 2; end
    if ~isfield(opts, 'kmax_auto'); opts.kmax_auto = 100; end
    if ~isfield(opts, 'hmin_adapt'); opts.hmin_adapt = 0; end
    if ~isfield(opts, 'maxit_adapt'); opts.maxit_adapt = 1e6; end

    L = b - a;
    fe_total = 0;

    usar_apriori = 0;
    M4 = NaN;
    if isfield(opts, 'M4')
        if isscalar(opts.M4) && isfinite(opts.M4) && opts.M4 >= 0
            usar_apriori = 1;
            M4 = opts.M4;
        end
    end

    % 1) Intento apriori si hay M4 finito
    if usar_apriori
        hreq = ((180.0 * tol) / (L * max(M4, eps)))^(1/4);
        n = ceil(L / hreq);
        if n < 2; n = 2; end
        if mod(n, 2) ~= 0; n = n + 1; end
        h = L / n;

        [Sn, fe1] = simpson_once_local(f, a, b, n);
        [S2n, fe2] = simpson_once_local(f, a, b, 2 * n);
        fe_total = fe_total + fe1 + fe2;

        err_est = abs(S2n - Sn) / 15.0;

        if err_est <= tol
            I = S2n;
            out = base_out('apriori', h, 2 * n, err_est, M4, fe_total);
            out.detalle = struct();
            out.detalle.n_apriori = n;
            out.detalle.h_apriori = h;
            out.detalle.bound_apriori = (L / 180.0) * (h^4) * M4;
            out.detalle.verificado = 1;
            return
        end

        % No alcanzo tol: continuar con halving global partiendo de 2n
        opts.n0 = 2 * n;
    end

    % 2) Halving global uniforme
    [Ig, hg, ng, errg, feg] = halving_global(f, a, b, tol, opts.n0, opts.kmax_auto);
    fe_total = fe_total + feg;
    if isfinite(errg) && errg <= tol
        I = Ig;
        out = base_out('uniforme_auto', hg, ng, errg, M4, fe_total);
        out.detalle = struct();
        out.detalle.k_iter = log2(ng / max(2, make_par(opts.n0)));
        return
    end

    % 3) Adaptativo (fallback)
    [Ia, infoA] = adaptativo_local(f, a, b, tol, opts.hmin_adapt, opts.maxit_adapt);
    fe_total = fe_total + infoA.fevals;
    I = Ia;

    out = base_out('adaptativo', NaN, NaN, NaN, M4, fe_total);
    out.detalle = infoA;
end

function n2 = make_par(n)
% Ajusta n a par, n >= 2
    n2 = floor(n);
    if n2 < 2; n2 = 2; end
    if mod(n2, 2) ~= 0; n2 = n2 + 1; end
end

function out = base_out(metodo, h, n, err_est, M4, fe_total)
% Arma struct base de salida
    out = struct();
    out.metodo = metodo;
    out.h = h;
    out.n = n;
    out.err_est = err_est;
    out.M4_used = M4;
    out.fevals_total = fe_total;
    out.detalle = struct();
end

function [S, fe] = simpson_once_local(f, a, b, n)
% Simpson 1/3 compuesto en malla uniforme con n par; retorna S y #evals
    if ~isscalar(n) || n < 2 || mod(n, 2) ~= 0 || floor(n) ~= n
        error('Simpson 1/3 requiere n par y n >= 2.');
    end
    h = (b - a) / n;

    fa = feval(f, a);
    fb = feval(f, b);
    fe = 2;

    s4 = 0;
    i = 1;
    while i <= n - 1
        xi = a + i * h;
        s4 = s4 + feval(f, xi);
        i = i + 2;
    end
    fe = fe + ceil((n - 1) / 2);

    s2 = 0;
    i = 2;
    while i <= n - 2
        xi = a + i * h;
        s2 = s2 + feval(f, xi);
        i = i + 2;
    end
    fe = fe + max(0, floor((n - 2) / 2));

    S = (h / 3.0) * (fa + fb + 4.0 * s4 + 2.0 * s2);
end

function [I, h, n, err_est, fevals] = halving_global(f, a, b, tol, n0, kmax)
% Halving global uniforme con estimacion p=4; retorna I,h,n,err,evals
    n = floor(n0);
    if n < 2; n = 2; end
    if mod(n, 2) ~= 0; n = n + 1; end
    kmax = floor(kmax);
    if kmax < 1; kmax = 1; end

    fevals = 0;

    [Sn, fe1] = simpson_once_local(f, a, b, n);
    fevals = fevals + fe1;

    err_est = Inf;
    k = 0;

    while k < kmax
        n2 = 2 * n;
        [S2n, fe2] = simpson_once_local(f, a, b, n2);
        fevals = fevals + fe2;

        err_est = abs(S2n - Sn) / 15.0;
        if err_est <= tol
            I = S2n;
            h = (b - a) / n2;
            n = n2;
            return
        end

        Sn = S2n;
        n = n2;
        k = k + 1;
    end

    I = Sn;
    h = (b - a) / n;
    warning('Halving: no se alcanzo tol en kmax iteraciones. err_est=%.3e', err_est);
end

function [I, info] = adaptativo_local(f, a, b, tol, hmin, maxit)
% Simpson 1/3 adaptativo iterativo (test 1/15 con correccion local)
    if nargin < 5 || isempty(hmin); hmin = 0; end
    if nargin < 6 || isempty(maxit); maxit = 1e6; end

    m = 0.5 * (a + b);
    fa = feval(f, a);
    fm = feval(f, m);
    fb = feval(f, b);
    fevals = 3;

    Sab = (b - a) * (fa + 4.0 * fm + fb) / 6.0;

    stack = [a, m, b, fa, fm, fb, Sab, tol];
    I = 0;
    n_accept = 0;
    n_split = 0;
    accepted = zeros(0, 3);
    err_loc = zeros(0, 1);

    it = 0;
    while ~isempty(stack)
        it = it + 1;
        if it > maxit
            warning('Adaptativo: se alcanzo maxit.');
            break
        end

        row = stack(end, :);
        stack(end, :) = [];

        ai = row(1); mi = row(2); bi = row(3);
        fai = row(4); fmi = row(5); fbi = row(6);
        Sab_i = row(7); tol_i = row(8);

        ci = 0.5 * (ai + mi);
        di = 0.5 * (mi + bi);
        fci = feval(f, ci);
        fdi = feval(f, di);
        fevals = fevals + 2;

        Sleft  = (mi - ai) * (fai + 4.0 * fci + fmi) / 6.0;
        Sright = (bi - mi) * (fmi + 4.0 * fdi + fbi) / 6.0;

        err = abs(Sleft + Sright - Sab_i) / 15.0;

        if err <= tol_i
            Scorr = Sleft + Sright + (Sleft + Sright - Sab_i) / 15.0;
            I = I + Scorr;
            n_accept = n_accept + 1;
            accepted(end+1, :) = [ai, mi, bi];
            err_loc(end+1, 1) = err;
        else
            if (hmin > 0) && ((mi - ai) <= hmin || (bi - mi) <= hmin)
                I = I + (Sleft + Sright);
                n_accept = n_accept + 1;
                accepted(end+1, :) = [ai, mi, bi];
                err_loc(end+1, 1) = err;
            else
                stack(end+1, :) = [mi, di, bi, fmi, fdi, fbi, Sright, tol_i / 2.0];
                stack(end+1, :) = [ai, ci, mi, fai, fci, fmi, Sleft,  tol_i / 2.0];
                n_split = n_split + 1;
            end
        end
    end

    info = struct();
    info.fevals = fevals;
    info.n_accept = n_accept;
    info.n_split = n_split;
    info.hmin = hmin;
    info.accepted = accepted;
    info.err_loc = err_loc;
end

