function [I, h, n, err_est, fevals] = trapecio_f_auto_min(f, a, b, tol, n0, kmax)
% trapecio_f_auto_min  Trapecios compuesto buscando el menor n que cumple tol.
% Inputs:
%   f     handle a funcion escalar, ej: @(x) sin(x)
%   a,b   limites con b > a
%   tol   tolerancia absoluta objetivo
%   n0    (opcional) n inicial (entero >= 1). Default 1
%   kmax  (opcional) max duplicaciones en etapa global. Default 20
% Outputs:
%   I        aproximacion final usando el n minimo que cumple (se devuelve I en 2*n)
%   h        paso final h = (b - a) / n, donde n es el minimo hallado
%   n        menor entero >=1 tal que |T_{2n} - T_n|/3 <= tol
%   err_est  estimacion de error asociada a I_{2n}: |T_{2n} - T_n|/3
%   fevals   cantidad total de evaluaciones de f
%
% Metodo:
%   1) Duplicar n hasta que |T_{2n}-T_n|/3 <= tol (orden p=2).
%   2) Busqueda fina: probar n = last_fail_n+1, ..., last_pass_n para encontrar
%      el menor n que satisface la tolerancia. Para cada candidato n_try se evalua
%      T_ntry e T_{2ntry} y se verifica |T_{2ntry}-T_{ntry}|/3 <= tol.
%
% Notas:
%   - f debe ser finita en los nodos.
%   - Esta variante devuelve I en el nivel 2*n (el usado en el test).
%   - Para cota apriori teorica: |E| <= (b-a)/12 * h^2 * M2, con M2 >= max|f''|.

    if nargin < 4
        error('Se requieren f, a, b y tol.');
    end
    if ~isa(f, 'function_handle')
        error('f debe ser function handle.');
    end
    if ~isscalar(a) || ~isscalar(b) || ~(b > a)
        error('Se requiere a < b.');
    end
    if ~isscalar(tol) || ~(tol > 0) || ~isfinite(tol)
        error('tol debe ser escalar positivo y finito.');
    end
    if nargin < 5 || isempty(n0)
        n0 = 1;
    end
    if nargin < 6 || isempty(kmax)
        kmax = 20;
    end

    n0 = floor(n0);
    if n0 < 1
        n0 = 1;
    end
    kmax = floor(kmax);
    if kmax < 1
        kmax = 1;
    end

    fevals = 0;

    % Etapa 1: duplicar n hasta cumplir tolerancia
    n = n0;
    [Tn, fe] = trapecio_once_(f, a, b, n);
    fevals = fevals + fe;

    passed = 0;
    k = 0;
    last_fail_n = 0;
    last_pass_n = 0;
    T2n_pass = NaN;
    err_at_pair = NaN;

    while k < kmax
        n2 = 2 * n;
        [T2n, fe2] = trapecio_once_(f, a, b, n2);
        fevals = fevals + fe2;

        err = abs(T2n - Tn) / 3.0;  % p=2 -> 2^p - 1 = 3

        if err <= tol
            passed = 1;
            last_fail_n = n;
            last_pass_n = n2;
            T2n_pass = T2n;
            err_at_pair = err;
            break
        end

        Tn = T2n;
        n = n2;
        k = k + 1;
    end

    if passed == 0
        % No alcanzo la tolerancia en kmax duplicaciones
        I = Tn;
        h = (b - a) / n;
        n = n;
        err_est = NaN;
        warning('trapecio_f_auto_min: no se alcanzo la tolerancia en kmax duplicaciones.');
        return
    end

    % Etapa 2: busqueda fina del menor n en (last_fail_n, last_pass_n]
    best_n = last_pass_n;
    best_I2 = T2n_pass;
    best_err = err_at_pair;

    n_try = last_fail_n + 1;
    while n_try <= last_pass_n
        [T_try, fea] = trapecio_once_(f, a, b, n_try);
        [T2_try, feb] = trapecio_once_(f, a, b, 2 * n_try);
        fevals = fevals + fea + feb;

        err_try = abs(T2_try - T_try) / 3.0;

        if err_try <= tol
            best_n = n_try;
            best_I2 = T2_try;
            best_err = err_try;
            break
        end

        n_try = n_try + 1;
    end

    n = best_n;
    h = (b - a) / n;
    I = best_I2;
    err_est = best_err;
end

function [T, fe] = trapecio_once_(f, a, b, n)
% trapecio_once_  Trapecios compuesto en [a,b] con n subintervalos.
% Devuelve T y el conteo aproximado de evaluaciones (n+1 nodos).

    if ~isscalar(n) || n < 1 || floor(n) ~= n
        error('n debe ser entero >= 1.');
    end

    h = (b - a) / n;

    fa = feval(f, a);
    fb = feval(f, b);
    if ~isfinite(fa) || ~isfinite(fb)
        error('Valores no finitos en los extremos.');
    end

    s = 0;
    i = 1;
    while i <= n - 1
        xi = a + i * h;
        fi = feval(f, xi);
        if ~isfinite(fi)
            error('Valor no finito de f en un nodo interior.');
        end
        s = s + fi;
        i = i + 1;
    end

    T = h * (0.5 * (fa + fb) + s);
    fe = n + 1;
end

