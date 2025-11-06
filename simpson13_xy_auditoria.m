function [I, info] = simpson13_xy_auditoria(x, y, a, b)
% simpson13_xy_auditoria  Integral en [a,b] usando Simpson 1/3 compuesto,
% con auditoria detallada y manejo de casos de borde.
%
% Inputs:
%   x  nodos estrictamente crecientes (malla uniformemente espaciada)
%   y  valores f(x) del mismo tamano que x
%   a  extremo izquierdo (x(1) <= a < b <= x(end))
%   b  extremo derecho
%
% Outputs:
%   I     aproximacion de la integral en [a,b]
%   info  struct con auditoria:
%         .a, .b
%         .h                      paso medio de la malla
%         .uniforme               1 si malla uniforme (dentro de tolerancia)
%         .interp_a, .interp_b    1 si se interpolo en el extremo
%         .left_idx, .right_idx   indices de x usados como interior
%         .X_interior, .Y_interior nodos y valores usados interiormente
%         .n_sub_interior         cantidad de subintervalos interiores
%         .blocks                 arreglo de bloques usados (ver abajo)
%         .n_blocks_trap          cantidad de bloques trapecio
%         .n_blocks_simp13        cantidad de bloques Simpson 1/3
%         .n_blocks_simp38        cantidad de bloques Simpson 3/8
%         .tol_x, .tol_h          tolerancias usadas
%
% info.blocks(k) tiene:
%   .tipo       'trapecio' | 'simpson13' | 'simpson38'
%   .x          vector de nodos del bloque
%   .y          vector de valores del bloque
%   .aporte     contribucion numerica de ese bloque
%   .indices    indices de x involucrados (si aplica)
%
% Notas/supuestos:
%   - x debe ser estrictamente creciente y aproximadamente uniforme.
%   - Simpson 1/3 requiere numero par de subintervalos en el tramo interior.
%   - Si a o b no son nodos, se aplica trapecio en el borde para alinear.
%   - Si el tramo interior queda con cantidad impar de subintervalos,
%     se usa 3/8 en los ultimos 3 subintervalos.
%   - Si el tramo interior queda con un solo subintervalo, se usa trapecio.
%
% Ejemplo de uso minimo:
%   x = 0:6:96;
%   y = [37.2 40.2 44.4 46.8 44.1 39.9 36.3 32.7 29.7 25.5 23.4 26.7 31.2 34.8 36.9 38.7 39.6];
%   [I, info] = simpson13_xy_auditoria(x, y, 0, 90);

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

    % Verificacion de uniformidad
    hx = diff(x);
    tol_h = 1e-10 * max(1, max(abs(hx)));
    if max(hx) - min(hx) > tol_h
        error('La malla x no es uniformemente espaciada dentro de la tolerancia.');
    end
    h = mean(hx);

    % Tolerancia para comparar a y b con nodos
    scaleX = max(1, max(abs(x)));
    tol_x = 1e-12 * scaleX;

    % Encontrar indices cercanos a a y b
    ia = find(x <= a + tol_x, 1, 'last');
    ib = find(x >= b - tol_x, 1, 'first');

    % Preparar auditoria
    blocks = struct('tipo', {}, 'x', {}, 'y', {}, 'aporte', {}, 'indices', {});
    n_trap = 0;
    n_s13 = 0;
    n_s38 = 0;

    % Caso: a y b en el mismo subintervalo -> un trapecio con interpolacion
    a_is_node = abs(a - x(ia)) <= tol_x;
    b_is_node = abs(b - x(ib)) <= tol_x;
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

        blk.tipo = 'trapecio';
        blk.x = [a; b];
        blk.y = [ya; yb];
        blk.aporte = I;
        blk.indices = [ia; ib];
        blocks(1) = blk;
        n_trap = 1;

        info = build_info(a, b, h, 1, ~a_is_node, ~b_is_node, ia, ib, [], [], 0, blocks, n_trap, n_s13, n_s38, tol_x, tol_h);
        return
    end

    I = 0;

    % Alinear extremo izquierdo con trapecio si a no es nodo
    if a_is_node
        left_idx = ia;
        interp_a = 0;
    else
        xa1 = x(ia);
        xa2 = x(ia+1);
        ya1 = y(ia);
        ya2 = y(ia+1);
        ya = ya1 + (ya2 - ya1) * (a - xa1) / (xa2 - xa1);
        aporte = 0.5 * (xa2 - a) * (ya + ya2);
        I = I + aporte;

        blk.tipo = 'trapecio';
        blk.x = [a; xa2];
        blk.y = [ya; ya2];
        blk.aporte = aporte;
        blk.indices = [ia; ia+1];
        blocks(end+1) = blk; %#ok<AGROW>
        n_trap = n_trap + 1;

        left_idx = ia + 1;
        interp_a = 1;
    end

    % Alinear extremo derecho con trapecio si b no es nodo
    if b_is_node
        right_idx = ib;
        interp_b = 0;
    else
        xb1 = x(ib-1);
        xb2 = x(ib);
        yb1 = y(ib-1);
        yb2 = y(ib);
        yb = yb1 + (yb2 - yb1) * (b - xb1) / (xb2 - xb1);
        aporte = 0.5 * (b - xb1) * (yb + yb1);
        I = I + aporte;

        blk.tipo = 'trapecio';
        blk.x = [xb1; b];
        blk.y = [yb1; yb];
        blk.aporte = aporte;
        blk.indices = [ib-1; ib];
        blocks(end+1) = blk; %#ok<AGROW>
        n_trap = n_trap + 1;

        right_idx = ib - 1;
        interp_b = 1;
    end

    % Si no quedan subintervalos interiores, terminar
    if right_idx <= left_idx
        info = build_info(a, b, h, 1, interp_a, interp_b, left_idx, right_idx, [], [], 0, blocks, n_trap, n_s13, n_s38, tol_x, tol_h);
        return
    end

    % Tramo interior
    Xint = x(left_idx:right_idx);
    Yint = y(left_idx:right_idx);
    nsub = right_idx - left_idx; % cantidad de subintervalos interiores

    % Caso especial: un solo subintervalo interior -> trapecio interior
    if nsub == 1
        aporte = 0.5 * (Xint(2) - Xint(1)) * (Yint(2) + Yint(1));
        I = I + aporte;

        blk.tipo = 'trapecio';
        blk.x = Xint(1:2);
        blk.y = Yint(1:2);
        blk.aporte = aporte;
        blk.indices = [left_idx; right_idx];
        blocks(end+1) = blk; %#ok<AGROW>
        n_trap = n_trap + 1;

        info = build_info(a, b, h, 1, interp_a, interp_b, left_idx, right_idx, Xint, Yint, nsub, blocks, n_trap, n_s13, n_s38, tol_x, tol_h);
        return
    end

    % Caso general
    if mod(nsub, 2) == 0
        % Par: Simpson 1/3 en todos los pares
        i = left_idx;
        while i + 2 <= right_idx
            aporte = (h/3) * (y(i) + 4*y(i+1) + y(i+2));
            I = I + aporte;

            blk.tipo = 'simpson13';
            blk.x = [x(i); x(i+1); x(i+2)];
            blk.y = [y(i); y(i+1); y(i+2)];
            blk.aporte = aporte;
            blk.indices = [i; i+1; i+2];
            blocks(end+1) = blk; %#ok<AGROW>
            n_s13 = n_s13 + 1;

            i = i + 2;
        end
    else
        % Impar: Simpson 1/3 en los primeros nsub-3, y 3/8 en los ultimos 3
        last_s13_end = left_idx + (nsub - 3);
        i = left_idx;
        while i + 2 <= last_s13_end
            aporte = (h/3) * (y(i) + 4*y(i+1) + y(i+2));
            I = I + aporte;

            blk.tipo = 'simpson13';
            blk.x = [x(i); x(i+1); x(i+2)];
            blk.y = [y(i); y(i+1); y(i+2)];
            blk.aporte = aporte;
            blk.indices = [i; i+1; i+2];
            blocks(end+1) = blk; %#ok<AGROW>
            n_s13 = n_s13 + 1;

            i = i + 2;
        end
        % Bloque final 3/8: cuatro nodos
        i0 = last_s13_end;
        i1 = i0 + 1;
        i2 = i0 + 2;
        i3 = i0 + 3;
        aporte = (3*h/8) * (y(i0) + 3*y(i1) + 3*y(i2) + y(i3));
        I = I + aporte;

        blk.tipo = 'simpson38';
        blk.x = [x(i0); x(i1); x(i2); x(i3)];
        blk.y = [y(i0); y(i1); y(i2); y(i3)];
        blk.aporte = aporte;
        blk.indices = [i0; i1; i2; i3];
        blocks(end+1) = blk; %#ok<AGROW>
        n_s38 = n_s38 + 1;
    end

    info = build_info(a, b, h, 1, interp_a, interp_b, left_idx, right_idx, Xint, Yint, nsub, blocks, n_trap, n_s13, n_s38, tol_x, tol_h);
end

function info = build_info(a, b, h, uniforme, interp_a, interp_b, left_idx, right_idx, Xint, Yint, nsub, blocks, n_trap, n_s13, n_s38, tol_x, tol_h)
% build_info  arma el struct de auditoria (funcion local)

    info = struct();
    info.a = a;
    info.b = b;
    info.h = h;
    info.uniforme = uniforme;
    info.interp_a = interp_a;
    info.interp_b = interp_b;
    info.left_idx = left_idx;
    info.right_idx = right_idx;
    info.X_interior = Xint;
    info.Y_interior = Yint;
    info.n_sub_interior = nsub;
    info.blocks = blocks;
    info.n_blocks_trap = n_trap;
    info.n_blocks_simp13 = n_s13;
    info.n_blocks_simp38 = n_s38;
    info.tol_x = tol_x;
    info.tol_h = tol_h;
end

