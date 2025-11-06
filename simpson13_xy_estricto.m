function I = simpson13_xy_estricto(x, y, a, b)
% simpson13_xy_estricto  Integral en [a,b] por Simpson 1/3 compuesto (estricto).
% Inputs:
%   x  nodos estrictamente crecientes y uniformemente espaciados
%   y  valores f(x) del mismo tamano que x
%   a  extremo izquierdo (debe coincidir con un nodo de x)
%   b  extremo derecho  (debe coincidir con un nodo de x)
% Output:
%   I  aproximacion de la integral en [a,b]
%
% Condiciones estrictas:
%   - x uniforme (misma separacion h dentro de tolerancia)
%   - a y b pertenecen a x (dentro de tolerancia)
%   - numero de subintervalos n = (indice(b)-indice(a)) es par y n >= 2
%
% Ejemplo minimo:
%   x = 0:6:36; y = sin(x*pi/180);
%   I = simpson13_xy_estricto(x, y, 0, 36);

    if nargin < 4
        error('Se requieren x, y, a, b.');
    end
    if ~isvector(x) || ~isvector(y) || numel(x) ~= numel(y)
        error('x e y deben ser vectores del mismo tamano.');
    end
    x = x(:);
    y = y(:);
    if any(~isfinite(x)) || any(~isfinite(y))
        error('x e y deben ser finitos.');
    end
    if any(diff(x) <= 0)
        error('x debe ser estrictamente creciente.');
    end

    hx = diff(x);
    tol_h = 1e-10 * max(1, max(abs(hx)));
    if max(hx) - min(hx) > tol_h
        error('x no es uniformemente espaciado dentro de la tolerancia.');
    end
    h = mean(hx);

    scaleX = max(1, max(abs(x)));
    tol_x = 1e-12 * scaleX;

    ia = find(abs(x - a) <= tol_x, 1, 'first');
    ib = find(abs(x - b) <= tol_x, 1, 'first');
    if isempty(ia) || isempty(ib)
        error('a y/o b no coinciden con nodos de x (modo estricto).');
    end
    if ~(ib > ia)
        error('Se requiere a < b y ambos en x.');
    end

    n = ib - ia; % numero de subintervalos
    if n < 2 || mod(n, 2) ~= 0
        error('Simpson 1/3 estricto requiere n par y n >= 2 en [a,b].');
    end

    yi = y(ia:ib);
    m = numel(yi);

    s4 = 0;
    k = 2;
    while k <= m - 1
        s4 = s4 + yi(k);
        k = k + 2;
    end

    s2 = 0;
    k = 3;
    while k <= m - 2
        s2 = s2 + yi(k);
        k = k + 2;
    end

    I = (h / 3) * (yi(1) + yi(m) + 4 * s4 + 2 * s2);
end

