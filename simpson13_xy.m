function I = simpson13_xy(x, y, a, b)
% simpson13_xy  Integral en [a,b] usando Simpson 1/3 compuesto con datos (x,y).
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
%   - Si a o b no son nodos, se usa trapecio en el extremo para alinear.
%   - Si los subintervalos interiores son impares, se usa 3/8 en los ultimos 3.
%   - Si solo queda un subintervalo interior, se usa trapecio en ese tramo.
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

    % Si a y b caen en el mismo subintervalo -> trapecio directo con interp.
    ia = find(x <= a + 1e-12*max(1,abs(a)), 1, 'last');
    ib = find(x >= b - 1e-12*max(1,abs(b)), 1, 'first');
    if ia == ib - 1
        ya = y(ia) + (y(ia+1) - y(ia)) * (a - x(ia)) / (x(ia+1) - x(ia));
        yb = y(ib-1) + (y(ib) - y(ib-1)) * (b - x(ib-1)) / (x(ib) - x(ib-1));
        I = 0.5 * (b - a) * (ya + yb);
        return
    end

    I = 0;

    % Alinear extremo izquierdo: si a no es nodo, trapecio [a, x(ia+1)]
    if abs(a - x(ia)) <= 1e-12*max(1,abs(a))
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

    % Alinear extremo derecho: si b no es nodo, trapecio [x(ib-1), b]
    if abs(b - x(ib)) <= 1e-12*max(1,abs(b))
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

    % Ahora integrar interiormente con Simpson 1/3 (y 3/8 si hiciera falta)
    if right_idx <= left_idx
        return
    end

    nsub = right_idx - left_idx; % numero de subintervalos interiores
    if nsub == 1
        % Un solo subintervalo interior -> trapecio
        I = I + 0.5 * (x(right_idx) - x(left_idx)) * (y(right_idx) + y(left_idx));
        return
    end

    % Caso general: nsub >= 2
    if mod(nsub, 2) == 0
        % Par: Simpson 1/3 en todo el interior
        yi = y(left_idx:right_idx);
        % pesos: 1,4,2,4,...,2,4,1
        s_par = yi(1) + yi(end);
        s4 = 0;
        s2 = 0;
        k = 2;
        while k <= numel(yi)-1
            s4 = s4 + yi(k);
            if k+1 <= numel(yi)-1
                s2 = s2 + yi(k+1);
            end
            k = k + 2;
        end
        I = I + (h/3) * (s_par + 4*s4 + 2*s2);
    else
        % Impar: Simpson 1/3 en los primeros nsub-3, y 3/8 en los ultimos 3
        if nsub >= 3
            % Parte Simpson 1/3
            cut_idx = left_idx + (nsub - 3);
            yi = y(left_idx:cut_idx);
            if numel(yi) >= 3
                s_par = yi(1) + yi(end);
                s4 = 0;
                s2 = 0;
                k = 2;
                while k <= numel(yi)-1
                    s4 = s4 + yi(k);
                    if k+1 <= numel(yi)-1
                        s2 = s2 + yi(k+1);
                    end
                    k = k + 2;
                end
                I = I + (h/3) * (s_par + 4*s4 + 2*s2);
            elseif numel(yi) == 2
                I = I + 0.5 * (x(cut_idx) - x(left_idx)) * (y(cut_idx) + y(left_idx));
            end
            % Bloque final 3/8 sobre 4 nodos: cut_idx...cut_idx+3
            i0 = cut_idx;
            i1 = cut_idx + 1;
            i2 = cut_idx + 2;
            i3 = cut_idx + 3;
            I = I + (3*h/8) * (y(i0) + 3*y(i1) + 3*y(i2) + y(i3));
        else
            % nsub == 1 (ya tratado) o 2 -> usar Simpson 1/3 si 2, trapecio si 1
            if nsub == 2
                yi = y(left_idx:right_idx);
                s_par = yi(1) + yi(end);
                s4 = yi(2);
                I = I + (h/3) * (s_par + 4*s4);
            else
                I = I + 0.5 * (x(right_idx) - x(left_idx)) * (y(right_idx) + y(left_idx));
            end
        end
    end
end

