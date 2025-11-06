function [I, h] = simpson13_f(f, a, b, n)
% simpson13_f  Integral de f en [a,b] por Simpson 1/3 compuesto.
% Inputs:
%   f  handle a funcion escalar, por ejemplo: @(x) sin(x)
%   a  extremo izquierdo del intervalo
%   b  extremo derecho del intervalo (b > a)
%   n  numero de subintervalos (entero par, n >= 2)
% Outputs:
%   I  aproximacion de la integral de f en [a,b]
%   h  paso de malla uniforme (h = (b - a) / n)
%
% Notas/supuestos:
%   - Requiere malla uniforme con n par (Simpson 1/3).
%   - f evaluable en [a,b] y sin singularidades en el intervalo.
%   - Orden global esperado O(h^4) si f es suficientemente suave.

    if nargin < 4
        error('Se requieren f, a, b y n.');
    end
    if ~isa(f, 'function_handle')
        error('f debe ser un function handle, por ejemplo @(x) x.^2.');
    end
    if ~isscalar(a) || ~isscalar(b) || ~(b > a)
        error('Se requiere a y b escalares con b > a.');
    end
    if ~isscalar(n) || n < 2 || mod(n, 2) ~= 0 || floor(n) ~= n
        error('n debe ser entero par y n >= 2.');
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

    I = (h / 3) * (fa + fb + 4 * s4 + 2 * s2);
end

