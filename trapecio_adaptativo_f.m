function [I, info] = trapecio_adaptativo_f(f, a, b, tol, hmin, maxit)
% trapecio_adaptativo_f  Trapecios adaptativo en [a,b] con control local 1/3.
% Inputs:
%   f     handle a funcion escalar, ej: @(x) sin(x)
%   a,b   limites con b > a
%   tol   tolerancia absoluta objetivo (por intervalo padre)
%   hmin  (opcional) ancho minimo permitido de subintervalo (default 0)
%   maxit (opcional) max iteraciones/subdivisiones (default 1e6)
% Outputs:
%   I     aproximacion de la integral en [a,b]
%   info  struct de auditoria con campos:
%         .fevals   evaluaciones totales de f
%         .n_accept cantidad de intervalos padre aceptados
%         .n_split  cantidad de subdivisiones realizadas
%         .hmin     hmin usado
%         .accepted Kx3 con filas [ai, mi, bi] de padres aceptados
%         .err_loc  Kx1 con err_i = |Tleft+Tright-Tab|/3 por intervalo
%
% Metodo:
%   - En el padre [a,b], trapecio T_ab = (b-a)*(f(a)+f(b))/2.
%   - Al partir m=(a+b)/2, hijos [a,m] y [m,b], con T_left y T_right.
%   - Estimador local p=2: err = |T_left + T_right - T_ab| / 3.
%   - Si err <= tol_local, se acepta con correccion de Richardson local:
%       T_corr = T_left + T_right + (T_left + T_right - T_ab) / 3
%     (mejora de orden, análogo al primer paso de Romberg).
%   - Si no, se subdivide repartiendo tol_local/2 a cada hijo.
%
% Advertencias:
%   - f debe ser finita en todos los nodos. Si hmin > 0 y no se puede
%     subdividir mas, se acepta T_left + T_right sin correccion.

    if nargin < 4
        error('Se requieren f, a, b, tol.');
    end
    if ~isa(f, 'function_handle')
        error('f debe ser function handle.');
    end
    if ~isscalar(a) || ~isscalar(b) || ~(b > a)
        error('Se requiere a < b.');
    end
    if ~isscalar(tol) || ~(tol > 0) || ~isfinite(tol)
        error('tol debe ser escalar positivo y finito.');
    end
    if nargin < 5 || isempty(hmin)
        hmin = 0;
    end
    if nargin < 6 || isempty(maxit)
        maxit = 1e6;
    end

    % Nodos iniciales y trapecio del padre
    m = 0.5 * (a + b);
    fa = feval(f, a);
    fm = feval(f, m);
    fb = feval(f, b);
    if ~isfinite(fa) || ~isfinite(fm) || ~isfinite(fb)
        error('Valores no finitos de f en nodos iniciales.');
    end
    fevals = 3;

    Tab = (b - a) * (fa + fb) * 0.5;

    % Pila con filas: [ai, mi, bi, fai, fmi, fbi, Tab_i, tol_i]
    stack = [a, m, b, fa, fm, fb, Tab, tol];

    I = 0;
    n_accept = 0;
    n_split = 0;
    accepted = zeros(0, 3);
    err_loc = zeros(0, 1);

    it = 0;
    while ~isempty(stack)
        it = it + 1;
        if it > maxit
            warning('Se alcanzo maxit; finalizando con lo acumulado.');
            break
        end

        row = stack(end, :);
        stack(end, :) = [];

        ai  = row(1);
        mi  = row(2);
        bi  = row(3);
        fai = row(4);
        fmi = row(5);
        fbi = row(6);
        Tab_i = row(7);
        tol_i = row(8);

        % Hijos y sus trapecios
        Tleft  = (mi - ai) * (fai + fmi) * 0.5;
        Tright = (bi - mi) * (fmi + fbi) * 0.5;

        err = abs(Tleft + Tright - Tab_i) / 3.0;

        if err <= tol_i
            % Aceptar con correccion local (Richardson)
            Tcorr = Tleft + Tright + (Tleft + Tright - Tab_i) / 3.0;
            I = I + Tcorr;
            n_accept = n_accept + 1;
            accepted(end+1, :) = [ai, mi, bi];
            err_loc(end+1, 1) = err;
        else
            % Control de hmin
            if (hmin > 0) && ((mi - ai) <= hmin || (bi - mi) <= hmin)
                I = I + (Tleft + Tright);
                n_accept = n_accept + 1;
                accepted(end+1, :) = [ai, mi, bi];
                err_loc(end+1, 1) = err;
            else
                % Subdividir: nuevos puntos y valores
                ci = 0.5 * (ai + mi);
                di = 0.5 * (mi + bi);
                fci = feval(f, ci);
                fdi = feval(f, di);
                if ~isfinite(fci) || ~isfinite(fdi)
                    error('Valor no finito de f en subdivision.');
                end
                fevals = fevals + 2;

                % Trapecios padres para hijos (usados como Tab_i en su test)
                Tab_left  = Tleft;
                Tab_right = Tright;

                % Empujar hijo derecho [mi, di, bi]
                stack(end+1, :) = [mi, di, bi, fmi, fdi, fbi, Tab_right, tol_i * 0.5];
                % Empujar hijo izquierdo [ai, ci, mi]
                stack(end+1, :) = [ai, ci, mi, fai, fci, fmi, Tab_left,  tol_i * 0.5];

                n_split = n_split + 1;
            end
        end
    end

    info = struct();
    info.fevals = fevals;
    info.n_accept = n_accept;
    info.n_split = n_split;
    info.hmin = hmin;
    info.accepted = accepted;
    info.err_loc = err_loc;
end

