function [I, h, n, err_est, fevals] = trapecio_f_auto(f, a, b, tol, n0, kmax)
% trapecio_f_auto  Trapecios compuesto con ajuste automatico de n por tolerancia.
% Inputs:
%   f     handle a funcion escalar, ej: @(x) sin(x)
%   a,b   limites con b > a
%   tol   tolerancia absoluta objetivo (scalar > 0)
%   n0    (opcional) n inicial (entero >= 1). Default 1
%   kmax  (opcional) maximo numero de duplicaciones. Default 20
% Outputs:
%   I        aproximacion final (usa malla uniforme)
%   h        paso final de malla, h = (b - a) / n
%   n        numero final de subintervalos (entero)
%   err_est  estimacion a posteriori = |T_{2n} - T_n| / 3  (orden p=2)
%   fevals   cantidad total de evaluaciones de f realizadas
%
% Notas:
%   - Metodo de refinamiento global (halving): duplica n hasta cumplir tol.
%   - Requiere f finita en todos los nodos.
%   - Cota a priori si se requiere justificar h teorico: |E| <= (b-a)/12 * h^2 * M2,
%     con M2 >= max |f''(x)| en [a,b].

    if nargin < 4
        error('Se requieren f, a, b, tol.');
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

    n = n0;
    [Tn, fe1] = trapecio_once_(f, a, b, n);
    fevals = fevals + fe1;

    err_est = Inf;
    k = 0;

    while k < kmax
        n2 = 2 * n;
        [T2n, fe2] = trapecio_once_(f, a, b, n2);
        fevals = fevals + fe2;

        err_est = abs(T2n - Tn) / 3.0;  % p = 2 -> divisor 2^p - 1 = 3

        if err_est <= tol
            I = T2n;
            n = n2;
            h = (b - a) / n;
            return
        end

        Tn = T2n;
        n = n2;
        k = k + 1;
    end

    % Si no alcanzo la tolerancia, devolver el ultimo estado y advertir
    I = Tn;
    h = (b - a) / n;
    warning('trapecio_f_auto: no se alcanzo la tolerancia en kmax duplicaciones. err_est=%.3e', err_est);
end

function [T, fe] = trapecio_once_(f, a, b, n)
% trapecio_once_  Trapecios compuesto en [a,b] con n subintervalos (malla uniforme).
% Devuelve T y el conteo simple de evaluaciones (n+1 nodos).

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

