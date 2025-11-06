function [I, info] = simpson38_adaptativo_f(f, a, b, tol, hmin, maxit)
% simpson38_adaptativo_f  Simpson 3/8 adaptativo en [a,b] con control 1/15.
% Inputs:
%   f     handle a funcion escalar, ej: @(x) sin(x)
%   a,b   limites con b > a
%   tol   tolerancia absoluta objetivo (por intervalo padre)
%   hmin  (opcional) ancho minimo permitido de subintervalo (default 0)
%   maxit (opcional) maximo de iteraciones/subdivisiones (default 1e6)
% Outputs:
%   I     aproximacion de la integral
%   info  struct de auditoria:
%         .fevals   cantidad total de evaluaciones de f
%         .n_accept cantidad de intervalos padre aceptados
%         .n_split  cantidad de subdivisiones realizadas
%         .hmin     hmin usado
%         .accepted Kx3 con [ai, mi, bi] de padres aceptados
%         .err_loc  Kx1 con err_i = |Sleft+Sright-Sab|/15 de cada padre aceptado
%
% Notas:
%   - En 3/8 se usan nodos a, a+(b-a)/3, a+2(b-a)/3, b.
%   - Al subdividir [a,b] en m=(a+b)/2, se reutilizan los nodos a+(b-a)/3 y
%     a+2(b-a)/3 en los hijos. Se evalua ademas en m, a+(b-a)/6 y b-(b-a)/6.
%   - Test de aceptacion: |Sleft+Sright-Sab|/15 <= tol_local. Si pasa, se
%     suma la correccion local de Richardson: + (Sleft+Sright-Sab)/15.
%   - Si hmin > 0 y no se puede subdividir mas, se acepta Sleft+Sright.

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

    % Nodos y 3/8 inicial en [a,b]
    h = (b - a) / 3;
    c1 = a + h;
    c2 = a + 2 * h;
    fa = feval(f, a);
    fc1 = feval(f, c1);
    fc2 = feval(f, c2);
    fb = feval(f, b);
    fevals = 4;
    Sab = (3 * h / 8) * (fa + 3 * fc1 + 3 * fc2 + fb);

    % Pila: [a, b, fa, fc1, fc2, fb, Sab, tol_local]
    stack = [a, b, fa, fc1, fc2, fb, Sab, tol];

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
        bi  = row(2);
        fai = row(3);
        fci1 = row(4); % f en ai + (bi-ai)/3
        fci2 = row(5); % f en ai + 2*(bi-ai)/3
        fbi = row(6);
        Sab_i = row(7);
        tol_i = row(8);

        % Geometria local
        hi = (bi - ai) / 3;
        m  = 0.5 * (ai + bi);

        % Nodos nuevos requeridos para 3/8 en cada mitad
        % Izquierda [ai, m]: nodos a, a+(m-a)/3, a+2(m-a)/3, m
        % Derecha   [m, bi]: nodos m, m+(b-m)/3, m+2(b-m)/3, b
        xL1 = ai + (m - ai) / 3;      % = ai + (bi - ai)/6
        xL2 = ai + 2 * (m - ai) / 3;  % = ai + (bi - ai)/3  (== c1)
        xR1 = m  + (bi - m) / 3;      % = ai + 2*(bi - ai)/3 (== c2)
        xR2 = m  + 2 * (bi - m) / 3;  % = bi - (bi - ai)/6

        fm  = feval(f, m);
        fL1 = feval(f, xL1);
        fR2 = feval(f, xR2);
        fevals = fevals + 3;

        % Sleft y Sright por 3/8 en cada mitad (reutilizando fci1, fci2)
        Sleft  = (3 * (m - ai) / 8) * (fai + 3 * fL1 + 3 * fci1 + fm);
        Sright = (3 * (bi - m) / 8) * (fm  + 3 * fci2 + 3 * fR2 + fbi);

        err = abs(Sleft + Sright - Sab_i) / 15;

        if err <= tol_i
            % Aceptar con correccion local
            Scorr = Sleft + Sright + (Sleft + Sright - Sab_i) / 15;
            I = I + Scorr;
            n_accept = n_accept + 1;
            accepted(end+1, :) = [ai, m, bi];
            err_loc(end+1, 1) = err;
        else
            % Control de hmin
            if (hmin > 0) && ((m - ai) <= hmin || (bi - m) <= hmin)
                I = I + (Sleft + Sright);
                n_accept = n_accept + 1;
                accepted(end+1, :) = [ai, m, bi];
                err_loc(end+1, 1) = err;
            else
                % Subdividir: push derecha y push izquierda repartiendose tol/2
                % Hijo derecho [m, bi]: extremos fm, fbi; nodos internos c1R=c2 (ya), c2R=xR2
                Sab_right = Sright;
                stack(end+1, :) = [m, bi, fm, fci2, fR2, fbi, Sab_right, tol_i / 2];
                % Hijo izquierdo [ai, m]: extremos fai, fm; nodos internos c1L=xL1, c2L=c1 (ya)
                Sab_left  = Sleft;
                stack(end+1, :) = [ai, m, fai, fL1, fci1, fm, Sab_left,  tol_i / 2];
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

