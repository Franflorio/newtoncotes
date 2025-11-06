function [I, h, n, err_est, fevals] = simpson38_f_auto(f, a, b, tol, n0, kmax)
% simpson38_f_auto  Simpson 3/8 compuesto con ajuste automatico por tolerancia.
% Inputs:
%   f     handle a funcion escalar, ej: @(x) sin(x)
%   a,b   limites con b > a
%   tol   tolerancia absoluta objetivo
%   n0    (opcional) n inicial, multiplo de 3 (default 12; se ajusta a multiplo de 3)
%   kmax  (opcional) maximo numero de refinamientos (default 20)
% Outputs:
%   I        aproximacion final de la integral
%   h        paso final de malla (h = (b - a)/n)
%   n        numero final de subintervalos (multiplo de 3)
%   err_est  estimacion de error |I_{2n}-I_n|/15
%   fevals   evaluaciones totales de f
%
% Notas:
%   - Mantiene n multiplo de 3 en todo momento.
%   - Estimacion de error por Richardson (orden 4), sin extrapolar el valor.
%   - f debe ser finita en los nodos.

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
        n0 = 12;
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

    n = n0;
    I_n = simpson38_once_(f, a, b, n);
    fevals = fevals + (n + 1);

    err_est = Inf;
    k = 0;

    while k < kmax
        n2 = 2 * n;  % sigue siendo multiplo de 3
        I_n2 = simpson38_once_(f, a, b, n2);
        fevals = fevals + (n2 + 1);

        err_est = abs(I_n2 - I_n) / 15;

        if ~isfinite(err_est)
            error('Estimacion de error no finita (posible f no finita en algun nodo).');
        end

        if err_est <= tol
            I = I_n2;
            n = n2;
            h = (b - a) / n;
            return
        end

        I_n = I_n2;
        n = n2;
        k = k + 1;
    end

    I = I_n;
    h = (b - a) / n;
    warning('No se alcanzo la tolerancia dentro de kmax refinamientos. err_est=%.3e', err_est);
end

function S = simpson38_once_(f, a, b, n)
% simpson38_once_  Simpson 3/8 compuesto en malla uniforme con n multiplo de 3.

    if ~isscalar(n) || n < 3 || mod(n, 3) ~= 0 || floor(n) ~= n
        error('Simpson 3/8 requiere n multiplo de 3 y n >= 3.');
    end

    h = (b - a) / n;

    fa = feval(f, a);
    fb = feval(f, b);
    if ~isfinite(fa) || ~isfinite(fb)
        error('Valores no finitos de f en los extremos.');
    end

    s2 = 0;  % indices i con mod(i,3)==0, internos
    s3 = 0;  % indices i con mod(i,3)~=0
    i = 1;
    while i <= n - 1
        xi = a + i * h;
        fi = feval(f, xi);
        if ~isfinite(fi)
            error('Valor no finito de f en un nodo interior.');
        end
        if mod(i, 3) == 0
            s2 = s2 + fi;
        else
            s3 = s3 + fi;
        end
        i = i + 1;
    end

    S = (3 * h / 8) * (fa + fb + 2 * s2 + 3 * s3);
end

