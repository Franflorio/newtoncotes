function I = trapecios_xy(x, y, a, b)
% trapecios_xy  Integral por trapecios en [a,b] usando x,y.
% Inputs:
%   x  vector de nodos (estrictamente creciente; no necesariamente uniforme)
%   y  vector de valores f(x) del mismo tamano que x
%   a  extremo izquierdo del intervalo de integracion (x(1) <= a < b <= x(end))
%   b  extremo derecho del intervalo de integracion
% Output:
%   I  aproximacion de la integral en [a,b] por regla de trapecios
%
% Notas/supuestos:
%   - x debe ser estrictamente creciente y del mismo tamano que y.
%   - Si a o b caen dentro de un subintervalo, se interpola linealmente.
%
% Alternativa:
%   -trapecios_vec(vector);
%   -[I, info] = trapecios_xy_auditoria(x,y,a,b);

% ======= VALIDACIONES INICIALES ========================
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
    if ~(isfinite(a) && isfinite(b) && b > a)
        error('Se requiere a < b y ambos finitos.');
    end
    if a < x(1) || b > x(end)
        error('a y b deben estar dentro de [x(1), x(end)].');
    end
%========= Tolerancia para decidir si a o b coinciden con algun nodo ======
    tol = 1e-12;

% Ubicar a y b dentro de la malla x
% Con esto identifico en que subintervalo cae a y en cual cae b

    ia_left = find(x <= a + tol, 1, 'last');
    ib_right = find(x >= b - tol, 1, 'first');

    if abs(a - x(ia_left)) <= tol
        ya = y(ia_left);
        left_interior_start = ia_left;
    else
        xa1 = x(ia_left);
        xa2 = x(ia_left + 1);
        ya1 = y(ia_left);
        ya2 = y(ia_left + 1);
        t = (a - xa1) / (xa2 - xa1);
        ya = ya1 + t * (ya2 - ya1);
        left_interior_start = ia_left;
    end

    if abs(b - x(ib_right)) <= tol
        yb = y(ib_right);
        right_interior_end = ib_right;
    else
        xb1 = x(ib_right - 1);
        xb2 = x(ib_right);
        yb1 = y(ib_right - 1);
        yb2 = y(ib_right);
        t = (b - xb1) / (xb2 - xb1);
        yb = yb1 + t * (yb2 - yb1);
        right_interior_end = ib_right;
    end

    idx_start = left_interior_start + 1;
    idx_end = right_interior_end - 1;

    if idx_start <= idx_end
        Xseg = [a; x(idx_start:idx_end); b];
        Yseg = [ya; y(idx_start:idx_end); yb];
    else
        Xseg = [a; b];
        Yseg = [ya; yb];
    end

    dX = Xseg(2:end) - Xseg(1:end-1);
    Ysum = Yseg(2:end) + Yseg(1:end-1);
    I = 0.5 * sum(dX .* Ysum);
end

