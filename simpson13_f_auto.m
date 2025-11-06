function [I, h, n, err_est, fevals] = simpson13_f_auto(f, a, b, tol, n0, kmax)
% simpson13_f_auto  Simpson 1/3 compuesto con autoajuste de n por tolerancia.
% Inputs:
%   f     handle a funcion escalar, por ejemplo: @(x) sin(x)
%   a,b   limites de integracion con b > a
%   tol   tolerancia deseada para el error (norma absoluta)
%   n0    (opcional) n inicial, entero par (default 10; si es impar se ajusta a par)
%   kmax  (opcional) maximo numero de refinamientos por halving de h (default 20)
% Outputs:
%   I        aproximacion final de la integral
%   h        paso final de malla (h = (b - a)/n)
%   n        numero final de subintervalos (par)
%   err_est  estimacion del error |I_{2n} - I_n| / 15
%   fevals   numero total de evaluaciones de f realizadas
%
% Notas/supuestos:
%   - Simpson 1/3 requiere n par; el algoritmo mantiene n par en todo momento.
%   - Control de error a posteriori con mallas anidadas (Richardson p=4).
%   - f debe ser evaluable y finita en [a,b]. Si hay singularidades, el
%     control de error puede no ser fiable.
%
% Ejemplo minimo:
%   f = @(x) sin(x);
%   [I,h,n,err,fev] = simpson13_f_auto(f, 0, pi, 1e-6)

    if nargin < 4
        error('Se requieren f, a, b y tol.');
    end
    if ~isa(f, 'function_handle')
        error('f debe ser un function handle.');
    end
    if ~isscalar(a) || ~isscalar(b) || ~(b > a)
        error('Se requiere a y b escalares con b > a.');
    end
    if ~isscalar(tol) || ~(tol > 0) || ~isfinite(tol)
        error('tol debe ser escalar positivo y finito.');
    end
    if nargin < 5
        n0 = 10;
    end
    if nargin < 6
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

    n = n0;
    I_n = simpson_once(f, a, b, n);
    fevals = fevals + (n + 1);

    k = 0;
    err_est = Inf;

    while k < kmax
        n2 = 2 * n;
        I_n2 = simpson_once(f, a, b, n2);
        fevals = fevals + (n2 + 1);

        err_est = abs(I_n2 - I_n) / 15;

        if ~isfinite(err_est)
            error('La estimacion de error no es finita (posibles valores no finitos de f).');
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

function S = simpson_once(f, a, b, n)
% simpson_once  Aplica Simpson 1/3 compuesto en una malla uniforme con n par.

    if mod(n, 2) ~= 0
        error('Simpson requiere n par.');
    end

    hloc = (b - a) / n;

    fa = feval(f, a);
    fb = feval(f, b);

    s4 = 0;
    i = 1;
    while i <= n - 1
        xi = a + i * hloc;
        s4 = s4 + feval(f, xi);
        i = i + 2;
    end

    s2 = 0;
    i = 2;
    while i <= n - 2
        xi = a + i * hloc;
        s2 = s2 + feval(f, xi);
        i = i + 2;
    end

    S = (hloc / 3) * (fa + fb + 4 * s4 + 2 * s2);
end

