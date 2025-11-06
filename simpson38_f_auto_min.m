function [I, h, n, err_est, fevals] = simpson38_f_auto_min(f, a, b, tol, n0, kmax)
% simpson38_f_auto_min
%   Simpson 3/8 compuesto con busqueda del menor n multiplo de 3 que cumple la tolerancia.
% Inputs:
%   f     handle a funcion escalar, ej: @(x) sin(x)
%   a,b   limites con b > a
%   tol   tolerancia absoluta objetivo
%   n0    (opcional) n inicial (multiplo de 3). Default 3; se ajusta a multiplo de 3
%   kmax  (opcional) maximo numero de duplicaciones globales. Default 20
% Outputs:
%   I        aproximacion final usando el n minimo que cumple (se devuelve I en 2*n)
%   h        paso final h = (b - a) / n, donde n es el menor multiplo de 3 hallado
%   n        menor numero de subintervalos multiplo de 3 tal que |I_{2n} - I_n|/15 <= tol
%   err_est  estimacion de error asociada a I_{2n}: err_est = |I_{2n} - I_n|/15
%   fevals   cantidad total de evaluaciones de f realizadas
%
% Metodo:
%   1) Duplicar n (refinamiento global) hasta que |I_{2n}-I_n|/15 <= tol.
%   2) Busqueda fina en pasos de 3 entre el ultimo n que no cumple y el primero que cumple,
%      eligiendo el menor n multiplo de 3 que satisface la tolerancia.
%
% Notas:
%   - Orden p=4, por eso el divisor 15 en la estimacion a posteriori.
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
        n0 = 3;
    end
    if nargin < 6 || isempty(kmax)
        kmax = 20;
    end

    n0 = floor(n0);
    if n0 < 3
        n0 = 3;
    end
    r = mod(n0, 3);
    if r ~= 0
        n0 = n0 + (3 - r);
    end
    kmax = floor(kmax);
    if kmax < 1
        kmax = 1;
    end

    fevals = 0;

    % Etapa 1: duplicar n hasta cumplir tolerancia
    n = n0;
    [Sn, fe1] = simpson38_once_(f, a, b, n);
    fevals = fevals + fe1;

    passed = 0;
    k = 0;
    last_fail_n = 0;
    last_pass_n = 0;
    S2n = NaN;

    while k < kmax
        n2 = 2 * n; % sigue siendo multiplo de 3
        [S2n_tmp, fe2] = simpson38_once_(f, a, b, n2);
        fevals = fevals + fe2;

        err = abs(S2n_tmp - Sn) / 15.0;

        if err <= tol
            passed = 1;
            last_fail_n = n;     % limite inferior de la busqueda fina
            last_pass_n = n2;    % primer valor que cumple
            S2n = S2n_tmp;
            break
        end

        Sn = S2n_tmp;
        n = n2;
        k = k + 1;
    end

    if passed == 0
        % No alcanzo la tolerancia dentro de kmax. Devolver ultimo estado.
        I = Sn;
        h = (b - a) / n;
        err_est = NaN;
        warning('No se alcanzo la tolerancia en kmax duplicaciones.');
        return
    end

    % Etapa 2: busqueda fina del menor n multiplo de 3 en (last_fail_n, last_pass_n]
    best_n = last_pass_n;
    best_I2 = S2n;
    best_err = abs(S2n - Sn) / 15.0;

    n_try = last_fail_n + 3;
    while n_try <= last_pass_n
        [S_try, fea] = simpson38_once_(f, a, b, n_try);
        [S2_try, feb] = simpson38_once_(f, a, b, 2 * n_try);
        fevals = fevals + fea + feb;

        err_try = abs(S2_try - S_try) / 15.0;

        if err_try <= tol
            best_n = n_try;
            best_I2 = S2_try;
            best_err = err_try;
            break
        end

        n_try = n_try + 3;
    end

    n = best_n;
    h = (b - a) / n;
    I = best_I2;
    err_est = best_err;
end

function [S, fe] = simpson38_once_(f, a, b, n)
% simpson38_once_  Simpson 3/8 compuesto en malla uniforme con n multiplo de 3.
% Devuelve S y el conteo simple de evaluaciones (n+1 nodos).
    if ~isscalar(n) || n < 3 || mod(n, 3) ~= 0 || floor(n) ~= n
        error('Simpson 3/8 requiere n multiplo de 3 y n >= 3.');
    end

    h = (b - a) / n;

    fa = feval(f, a);
    fb = feval(f, b);

    s2 = 0; % nodos internos i con mod(i,3)==0
    s3 = 0; % nodos internos i con mod(i,3)~=0

    i = 1;
    while i <= n - 1
        xi = a + i * h;
        fi = feval(f, xi);
        if mod(i, 3) == 0
            s2 = s2 + fi;
        else
            s3 = s3 + fi;
        end
        i = i + 1;
    end

    S = (3.0 * h / 8.0) * (fa + fb + 2.0 * s2 + 3.0 * s3);

    fe = n + 1;
end

