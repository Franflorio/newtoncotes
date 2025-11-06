function [I, h] = trapecio_f(f, a, b, n)
% trapecio_f  Regla de trapecios compuesta sobre [a,b] con n subintervalos.
% Inputs:
%   f  handle a funcion escalar, ej: @(x) sin(x)
%   a  limite inferior (escalar)
%   b  limite superior (escalar), b > a
%   n  cantidad de subintervalos (entero >= 1)
% Outputs:
%   I  aproximacion integral de f en [a,b] por trapecios compuesto
%   h  paso de malla, h = (b - a) / n
%
% Notas:
%   - Requiere f evaluable y finita en todos los nodos.
%   - Para cota apriori: si f'' es continua y |f''| <= M2,
%     entonces |E| <= (b-a)/12 * h^2 * M2.

    % Chequeos basicos
    if nargin < 4
        error('Se requieren f, a, b, n.');
    end
    if ~isa(f, 'function_handle')
        error('f debe ser un function handle.');
    end
    if ~isscalar(a) || ~isscalar(b) || ~(b > a)
        error('Se requiere a < b.');
    end
    if ~isscalar(n) || ~isfinite(n) || n < 1 || floor(n) ~= n
        error('n debe ser entero >= 1.');
    end

    h = (b - a) / n;

    % Extremos
    fa = feval(f, a);
    fb = feval(f, b);
    if ~isfinite(fa) || ~isfinite(fb)
        error('Valores no finitos en los extremos.');
    end

    % Suma de nodos internos
    s = 0;
    i = 1;
    while i <= n - 1
        xi = a + i * h;
        fi = feval(f, xi);
        if ~isfinite(fi)
            error('Valor no finito en un nodo interior.');
        end
        s = s + fi;
        i = i + 1;
    end

    % Formula compuesta
    I = h * (0.5 * (fa + fb) + s);
end

