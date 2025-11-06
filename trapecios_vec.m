function I = trapecios_vec(y)
% trapecios_vec  Integral por trapecios a partir de un vector y.
% Inputs:
%   y  vector con valores de f(x_i) en nodos equiespaciados (incluye a y b)
% Output:
%   I  aproximacion de la integral en [a,b]
%
% Notas/supuestos:
%   - y corresponde a n+1 nodos equiespaciados entre a y b.
%   - Este m-file pide a y b por teclado.
%
% Ejemplo rapido:
%   y = sin(linspace(0, pi, 101));
%   I = trapecios_vec(y);  % ingresar a=0, b=pi cuando lo pida
%   % valor exacto: 2

    if nargin < 1
        error('Se requiere el vector y.');
    end

    if ~isvector(y)
        error('y debe ser un vector.');
    end

    if any(~isfinite(y))
        error('y debe contener solo valores finitos.');
    end

    n = numel(y);
    if n < 2
        error('y debe tener al menos 2 elementos.');
    end

    a = input('Ingrese a (extremo izquierdo): ');
    b = input('Ingrese b (extremo derecho): ');

    if ~isscalar(a) || ~isscalar(b) || ~isfinite(a) || ~isfinite(b)
        error('a y b deben ser escalares finitos.');
    end

    if ~(b > a)
        error('Se requiere b > a.');
    end

    h = (b - a) / (n - 1);
    I = h * (sum(y) - 0.5 * (y(1) + y(end)));
end

