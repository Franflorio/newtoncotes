function [I, info] = trapecios_xy_auditoria(x, y, a, b)
% trapecios_xy_auditoria  Integral por trapecios en [a,b] con datos x,y.
% Inputs:
%   x  vector de nodos estrictamente creciente (no requiere malla uniforme)
%   y  vector con f(x) del mismo tamano que x
%   a  extremo izquierdo de integracion (x(1) <= a < b <= x(end))
%   b  extremo derecho de integracion
% Outputs:
%   I     aproximacion integral en [a,b] por regla de trapecios compuesta
%   info  struct con datos de auditoria:
%         .a, .b, .ya, .yb
%         .interp_a, .interp_b  (1 si interpolo en el extremo, 0 si uso nodo)
%         .idx_interior  indices de x usados como interiores
%         .Xseg, .Yseg   nodos usados incluyendo a y b
%         .h             anchos por subintervalo (diff(Xseg))
%         .sumas         contribucion de cada trapecio
%         .n_sub         cantidad de subintervalos
%         .uniforme      1 si todos los h son iguales dentro de tolerancia
%         .h_min, .h_max, .h_med
%         .tol           tolerancia numerica usada
%
% Notas/supuestos:
%   - x debe ser estrictamente creciente, y del mismo tamano que y.
%   - a y b deben estar dentro de [x(1), x(end)] y cumplir a < b.
%   - Si a o b no son nodos, se interpola linealmente en el extremo.
%   - No usa toolboxes. ASCII puro.

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

    scale = max(1, max(abs(x)));
    tol = 1e-12 * scale;

    ia_left = find(x <= a + tol, 1, 'last');
    ib_right = find(x >= b - tol, 1, 'first');

    if abs(a - x(ia_left)) <= tol
        ya = y(ia_left);
        interp_a = 0;
    else
        xa1 = x(ia_left);
        xa2 = x(ia_left + 1);
        ya1 = y(ia_left);
        ya2 = y(ia_left + 1);
        t = (a - xa1) / (xa2 - xa1);
        ya = ya1 + t * (ya2 - ya1);
        interp_a = 1;
    end

    if abs(b - x(ib_right)) <= tol
        yb = y(ib_right);
        interp_b = 0;
    else
        xb1 = x(ib_right - 1);
        xb2 = x(ib_right);
        yb1 = y(ib_right - 1);
        yb2 = y(ib_right);
        t = (b - xb1) / (xb2 - xb1);
        yb = yb1 + t * (yb2 - yb1);
        interp_b = 1;
    end

    idx_start = ia_left + 1;
    idx_end = ib_right - 1;

    if idx_start <= idx_end
        idx_interior = (idx_start:idx_end).';
        Xseg = [a; x(idx_interior); b];
        Yseg = [ya; y(idx_interior); yb];
    else
        idx_interior = [];
        Xseg = [a; b];
        Yseg = [ya; yb];
    end

    h = Xseg(2:end) - Xseg(1:end-1);
    sumas = 0.5 * h .* (Yseg(2:end) + Yseg(1:end-1));
    I = sum(sumas);

    if numel(h) >= 1
        tolH = 1e-12 * max(1, max(abs(h)));
        uniforme = (max(h) - min(h) <= tolH);
        h_min = min(h);
        h_max = max(h);
        h_med = mean(h);
    else
        uniforme = 1;
        h_min = b - a;
        h_max = b - a;
        h_med = b - a;
    end

    info = struct();
    info.a = a;
    info.b = b;
    info.ya = ya;
    info.yb = yb;
    info.interp_a = interp_a;
    info.interp_b = interp_b;
    info.idx_interior = idx_interior;
    info.Xseg = Xseg;
    info.Yseg = Yseg;
    info.h = h;
    info.sumas = sumas;
    info.n_sub = numel(h);
    info.uniforme = uniforme;
    info.h_min = h_min;
    info.h_max = h_max;
    info.h_med = h_med;
    info.tol = tol;
end

