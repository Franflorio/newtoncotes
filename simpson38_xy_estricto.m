function I = simpson38_xy_estricto(x, y, a, b)
% simpson38_xy_estricto  Integral en [a,b] por Simpson 3/8 compuesto (estricto).
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
%   - numero de subintervalos n = (indice(b)-indice(a)) multiplo de 3 y n >= 3
%
% Ejemplo minimo:
%   x = 0:6:36; y = sin(x*pi/180);
%   I = simpson38_xy_estricto(x, y, 0, 36);

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
    if n < 3 || mod(n, 3) ~= 0
        error('Simpson 3/8 estricto requiere n multiplo de 3 y n >= 3 en [a,b].');
    end

    I = 0;
    i = ia;
    while i + 3 <= ib
        I = I + (3 * h / 8) * (y(i) + 3 * y(i + 1) + 3 * y(i + 2) + y(i + 3));
        i = i + 3;
    end
end

