function [I, h, n, bound] = simpson13_f_apriori(f, a, b, tol, M4)
% simpson13_f_apriori  Simpson 1/3 compuesto con seleccion a priori de h y n.
% Inputs:
%   f    handle a funcion escalar, ej: @(x) sin(x)
%   a    extremo izquierdo
%   b    extremo derecho (b > a)
%   tol  tolerancia absoluta deseada para la cota teorica de error
%   M4   cota superior de |f''''(x)| en [a,b] (escalar >= 0)
% Outputs:
%   I      aproximacion de la integral de f en [a,b]
%   h      paso final de malla (h = (b - a)/n)
%   n      numero final de subintervalos (par)
%   bound  cota teorica del error con el h elegido: (b-a)/180 * h^4 * M4
%
% Notas/supuestos:
%   - Requiere M4 >= max|f''''| en [a,b]. Si M4 se subestima, no hay garantia.
%   - Si M4 == 0 (f polinomio de grado <= 3), Simpson es exacto: usa n = 2.
%   - f debe ser evaluable y finita en los nodos.
%   - No usa toolboxes. ASCII puro.
%
% Ejemplo rapido (ejercicio del fluido):
%   m = 10;
%   f = @(v) m ./ (v.^1.5);
%   a = 5; b = 19; tol = 1e-3;
%   M4 = (945*m/16) * a^(-11/2); % max en el extremo izquierdo para v^{-11/2}
%   [I,h,n,bound] = simpson13_f_apriori(f, a, b, tol, M4)

    if nargin < 5
        error('Se requieren f, a, b, tol y M4.');
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
    if ~isscalar(M4) || ~(M4 >= 0) || ~isfinite(M4)
        error('M4 debe ser una cota escalar no negativa y finita.');
    end

    L = b - a;

    if M4 == 0
        n = 2;  % Simpson exacto si f es de grado <= 3
        h = L / n;
    else
        hreq = (180 * tol / (L * M4))^(1/4);
        if ~(isfinite(hreq)) || hreq <= 0
            error('No se pudo calcular h requerido (verifique M4 y tol).');
        end
        n = ceil(L / hreq);
        if n < 2
            n = 2;
        end
        if mod(n, 2) ~= 0
            n = n + 1;
        end
        h = L / n;
    end

    % Simpson 1/3 compuesto
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

    I = (h / 3) * (fa + fb + 4 * s4 + 2 * s2);

    bound = (L / 180) * (h^4) * M4;
end

