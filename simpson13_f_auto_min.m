function [I, h, n, err_est, fevals] = simpson13_f_auto_min(f, a, b, tol, n0, kmax)
% simpson13_f_auto_min  Simpson 1/3 compuesto con busqueda del menor n par que cumple tol.
% Inputs:
%   f     handle a funcion escalar, ej: @(x) sin(x)
%   a,b   limites con b > a
%   tol   tolerancia absoluta objetivo
%   n0    (opcional) n inicial (par). Default 2; si es impar se ajusta a par
%   kmax  (opcional) maximo numero de duplicaciones globales. Default 20
% Outputs:
%   I        aproximacion final usando el n minimo que cumple (se devuelve I en 2*n)
%   h        paso final de malla h = (b - a) / n  (n es el menor par que cumple)
%   n        menor numero de subintervalos PAR tal que |I_{2n} - I_n|/15 <= tol
%   err_est  estimacion de error asociada a I_{2n}: err_est = |I_{2n} - I_n|/15
%   fevals   cantidad total de evaluaciones de f realizadas
%
% Metodo:
%   1) Halving global: duplicar n hasta que |I_{2n}-I_n|/15 <= tol.
%   2) Busqueda fina: entre el ultimo n fallido y el n que paso, probar n pares
%      crecientes y elegir el menor que cumple. Para cada candidato n_try se evalua
%      I_ntry e I_{2ntry} y se verifica |I_{2ntry}-I_ntry|/15 <= tol.
%
% Notas:
%   - f debe ser finita en todos los nodos de integracion.
%   - No usa toolboxes. ASCII puro. Una funcion por archivo.

    if nargin < 4
        error('Se requieren f, a, b y tol.');
    end
    if ~isa(f, 'function_handle')
        error('f debe ser un function handle.');
    end
    if ~isscalar(a) || ~isscalar(b) || ~(b > a)
        error('Se requiere a < b.');
    end
    if ~isscalar(tol) || ~(tol > 0) || ~isfinite(tol)
        error('tol debe ser escalar positivo y finito.');
    end
    if nargin < 5 || isempty(n0)
        n0 = 2;
    end
    if nargin < 6 || isempty(kmax)
        kmax = 20;
    end

    n0 = floor(n0);
    if n0 < 2
        n0 = 2;
    end
    if mod(n0, 2) ~= 0
        n0 = n0 + 1;
    end
    kmax = floor(kmax);
    if kmax < 1
        kmax = 1;
    end

    fevals = 0;

    % Etapa 1: duplicar n hasta cumplir tolerancia
    n = n0;
    [Sn, fe] = simpson_once_(f, a, b, n);
    fevals = fevals + fe;

    passed = 0;
    k = 0;
    last_fail_n = 0;
    last_pass_n = 0;
    S2n = NaN;

    while k < kmax
        n2 = 2 * n;
        [S2n_tmp, fe2] = simpson_once_(f, a, b, n2);
        fevals = fevals + fe2;

        err = abs(S2n_tmp - Sn) / 15.0;

        if err <= tol
            passed = 1;
            last_fail_n = n;    % este n es el limite inferior de la busqueda fina
            last_pass_n = n2;   % este 2n pasa
            S2n = S2n_tmp;
            break
        end

        % continuar duplicando
        Sn = S2n_tmp;
        n = n2;
        k = k + 1;
    end

    if passed == 0
        % No alcanzo la tolerancia dentro de kmax. Se devuelve el ultimo estado.
        I = Sn;
        h = (b - a) / n;
        err_est = NaN;
        warning('No se alcanzo la tolerancia en kmax duplicaciones.');
        return
    end

    % Etapa 2: busqueda fina del menor n par en (last_fail_n, last_pass_n]
    best_n = last_pass_n;
    best_I2 = S2n;
    best_err = abs(S2n - Sn) / 15.0;

    n_try = last_fail_n + 2;
    while n_try <= last_pass_n
        [S_try, fe1] = simpson_once_(f, a, b, n_try);
        [S2_try, fe2] = simpson_once_(f, a, b, 2 * n_try);
        fevals = fevals + fe1 + fe2;

        err_try = abs(S2_try - S_try) / 15.0;

        if err_try <= tol
            best_n = n_try;
            best_I2 = S2_try;
            best_err = err_try;
            break
        end

        n_try = n_try + 2;
    end

    n = best_n;
    h = (b - a) / n;
    I = best_I2;
    err_est = best_err;
end

function [S, fe] = simpson_once_(f, a, b, n)
% Simpson 1/3 compuesto en malla uniforme con n par. Devuelve S y #evals.
    if ~isscalar(n) || n < 2 || mod(n, 2) ~= 0 || floor(n) ~= n
        error('Simpson 1/3 requiere n par y n >= 2.');
    end

    h = (b - a) / n;

    fa = feval(f, a);
    fb = feval(f, b);

    s4 = 0;
    i = 1;
    while i <= n - 1
        xi = a + i * h;
        s4 = s4 + feval(f, xi);
        i = i + 2;
    end

    s2 = 0;
    i = 2;
    while i <= n - 2
        xi = a + i * h;
        s2 = s2 + feval(f, xi);
        i = i + 2;
    end

    S = (h / 3.0) * (fa + fb + 4.0 * s4 + 2.0 * s2);

    % conteo simple: n+1 evaluaciones nuevas
    fe = n + 1;
end

