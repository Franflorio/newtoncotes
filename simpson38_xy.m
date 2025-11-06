function I = simpson38_xy(x, y, a, b)
% simpson38_xy  Integral en [a,b] usando Simpson 3/8 compuesto con datos (x,y).
% Inputs:
%   x  nodos estrictamente crecientes (malla uniformemente espaciada)
%   y  valores f(x) del mismo tamano que x
%   a  extremo izquierdo (x(1) <= a < b <= x(end))
%   b  extremo derecho
% Output:
%   I  aproximacion de la integral en [a,b]
%
% Notas/supuestos:
%   - x debe ser estrictamente creciente y aproximadamente uniforme.
%   - El metodo 3/8 trabaja con bloques de 3 subintervalos (4 nodos).
%   - Si a o b no son nodos, se usa trapecio en el borde para alinear.
%   - Si el numero de subintervalos interiores no es multiplo de 3,
%     se combina con Simpson 1/3 (2 subintervalos) o trapecio (1).
%   - Requiere f con regularidad suficiente para la precision esperada.

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

    % Verificar uniformidad de x
    hx = diff(x);
    tolH = 1e-10 * max(1, max(abs(hx)));
    if max(hx) - min(hx) > tolH
        error('La malla x no es uniformemente espaciada dentro de la tolerancia.');
    end
    h = mean(hx);

    % Tolerancia para comparar a y b con nodos
    scaleX = max(1, max(abs(x)));
    tolX = 1e-12 * scaleX;

    ia = find(x <= a + tolX, 1, 'last');
    ib = find(x >= b - tolX, 1, 'first');

    a_is_node = abs(a - x(ia)) <= tolX;
    b_is_node = abs(b - x(ib)) <= tolX;

    % Caso: a y b en el mismo subintervalo -> trapecio con interpolacion
    if ia == ib - 1
        if a_is_node
            ya = y(ia);
        else
            ya = y(ia) + (y(ia+1) - y(ia)) * (a - x(ia)) / (x(ia+1) - x(ia));
        end
        if b_is_node
            yb = y(ib);
        else
            yb = y(ib-1) + (y(ib) - y(ib-1)) * (b - x(ib-1)) / (x(ib) - x(ib-1));
        end
        I = 0.5 * (b - a) * (ya + yb);
        return
    end

    I = 0;

    % Alinear extremo izquierdo con trapecio si a no es nodo
    if a_is_node
        left_idx = ia;
    else
        xa1 = x(ia);
        xa2 = x(ia+1);
        ya1 = y(ia);
        ya2 = y(ia+1);
        ya = ya1 + (ya2 - ya1) * (a - xa1) / (xa2 - xa1);
        I = I + 0.5 * (xa2 - a) * (ya + ya2);
        left_idx = ia + 1;
    end

    % Alinear extremo derecho con trapecio si b no es nodo
    if b_is_node
        right_idx = ib;
    else
        xb1 = x(ib-1);
        xb2 = x(ib);
        yb1 = y(ib-1);
        yb2 = y(ib);
        yb = yb1 + (yb2 - yb1) * (b - xb1) / (xb2 - xb1);
        I = I + 0.5 * (b - xb1) * (yb + yb1);
        right_idx = ib - 1;
    end

    % Si no quedan subintervalos interiores, terminar
    if right_idx <= left_idx
        return
    end

    % Tramo interior
    nsub = right_idx - left_idx; % cantidad de subintervalos interiores

    % Casos pequenos
    if nsub == 1
        I = I + 0.5 * (x(right_idx) - x(left_idx)) * (y(right_idx) + y(left_idx));
        return
    end
    if nsub == 2
        I = I + (h/3) * (y(left_idx) + 4*y(left_idx+1) + y(left_idx+2));
        return
    end

    % Caso general nsub >= 3
    r = mod(nsub, 3);

    if r == 0
        % Solo bloques 3/8
        i = left_idx;
        while i + 3 <= right_idx
            I = I + (3*h/8) * (y(i) + 3*y(i+1) + 3*y(i+2) + y(i+3));
            i = i + 3;
        end

    elseif r == 2
        % Bloques 3/8 hasta dejar 2 subintervalos al final, que se cubren con 1/3
        end38 = left_idx + (nsub - 2);
        i = left_idx;
        while i + 3 <= end38
            I = I + (3*h/8) * (y(i) + 3*y(i+1) + 3*y(i+2) + y(i+3));
            i = i + 3;
        end
        % Ultimos 2 subintervalos -> Simpson 1/3
        I = I + (h/3) * (y(end38) + 4*y(end38+1) + y(end38+2));

    else
        % r == 1
        % Bloques 3/8 hasta dejar 4 subintervalos al final,
        % que se cubren con dos Simpson 1/3 consecutivos
        if nsub < 4
            % Seguridad: si nsub=1 o 2 ya cubierto arriba; si nsub=3 usar 3/8 directo
            if nsub == 3
                I = I + (3*h/8) * (y(left_idx) + 3*y(left_idx+1) + 3*y(left_idx+2) + y(left_idx+3));
            else
                I = I + 0.5 * (x(right_idx) - x(left_idx)) * (y(right_idx) + y(left_idx));
            end
            return
        end
        end38 = left_idx + (nsub - 4);
        i = left_idx;
        while i + 3 <= end38
            I = I + (3*h/8) * (y(i) + 3*y(i+1) + 3*y(i+2) + y(i+3));
            i = i + 3;
        end
        % Ultimos 4 subintervalos -> dos Simpson 1/3 que comparten el nodo central
        % Primer 1/3: nodos i, i+1, i+2
        I = I + (h/3) * (y(i) + 4*y(i+1) + y(i+2));
        % Segundo 1/3: nodos i+2, i+3, i+4 (i+4 == right_idx)
        I = I + (h/3) * (y(i+2) + 4*y(i+3) + y(i+4));
    end
end

