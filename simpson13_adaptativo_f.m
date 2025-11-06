function [I, info] = simpson13_adaptativo_f(f, a, b, tol, hmin, maxit)
% simpson13_adaptativo_f  Simpson 1/3 adaptativo en [a,b] con control 1/15.
% Inputs:
%   f     handle a funcion escalar, ej: @(x) sin(x)
%   a,b   limites con b > a
%   tol   tolerancia absoluta objetivo para el error local compuesto
%   hmin  (opcional) minimo ancho de subintervalo permitido (default 0)
%   maxit (opcional) maximo numero de iteraciones/subdivisiones (default 1e6)
% Outputs:
%   I     aproximacion de la integral
%   info  struct con auditoria:
%         .fevals   cantidad de evaluaciones de f
%         .n_accept cantidad de intervalos aceptados
%         .n_split  cantidad de subdivisiones realizadas
%         .hmin     hmin usado
%         .accepted matriz Kx3 con [ai, mi, bi] aceptados
%         .err_loc  vector Kx1 con estimaciones locales err_i = |Sli+Sri-Sab|/15
%
% Metodo (segun apunte):
%   1) S(a,b) = (b-a)/6*(f(a)+4*f((a+b)/2)+f(b)).
%   2) Partir en m=(a+b)/2, computar S(a,m) y S(m,b) reutilizando f.
%   3) Test: |S(a,m)+S(m,b)-S(a,b)|/15 <= tol? Si: aceptar con correccion
%      de Richardson local: S_corr = S(a,m)+S(m,b) + (S(a,m)+S(m,b)-S(a,b))/15.
%      Si no: subdividir y repartir tol/2 a cada hijo.
%   4) Repetir hasta cumplir el test o alcanzar hmin/maxit.
%
% Notas:
%   - Implementacion iterativa con pila para cumplir "una funcion por archivo".
%   - Si hmin > 0 y ya no se puede subdividir, se acepta la suma Sl+Sr sin
%     correccion y se continua (se prioriza terminar).
%   - No usa toolboxes. ASCII puro.

    if nargin < 4
        error('Se requieren f, a, b, tol.');
    end
    if ~isa(f, 'function_handle')
        error('f debe ser un function handle.');
    end
    if ~isscalar(a) || ~isscalar(b) || ~(b > a)
        error('Se requiere a < b.');
    end
    if ~isscalar(tol) || ~(tol > 0) || ~isfinite(tol)
        error('tol debe ser escalar positivo finito.');
    end
    if nargin < 5 || isempty(hmin)
        hmin = 0;
    end
    if nargin < 6 || isempty(maxit)
        maxit = 1e6;
    end

    % Evaluaciones iniciales
    m = 0.5 * (a + b);
    fa = feval(f, a);
    fm = feval(f, m);
    fb = feval(f, b);
    fevals = 3;

    Sab = (b - a) * (fa + 4 * fm + fb) / 6;

    % Pila con filas: [a, m, b, fa, fm, fb, Sab, tol_local]
    stack = [a, m, b, fa, fm, fb, Sab, tol];

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

        ai = row(1);
        mi = row(2);
        bi = row(3);
        fai = row(4);
        fmi = row(5);
        fbi = row(6);
        Sab_i = row(7);
        tol_i = row(8);

        % Subdividir el intervalo actual en dos mitades
        ci = 0.5 * (ai + mi);
        di = 0.5 * (mi + bi);

        fci = feval(f, ci);
        fdi = feval(f, di);
        fevals = fevals + 2;

        Sleft  = (mi - ai) * (fai + 4 * fci + fmi) / 6;
        Sright = (bi - mi) * (fmi + 4 * fdi + fbi) / 6;

        err = abs(Sleft + Sright - Sab_i) / 15;

        if err <= tol_i
            % Aceptar con correccion local de Richardson
            Scorr = Sleft + Sright + (Sleft + Sright - Sab_i) / 15;
            I = I + Scorr;
            n_accept = n_accept + 1;
            accepted(end+1, :) = [ai, mi, bi];
            err_loc(end+1, 1) = err;
        else
            % Verificar hmin
            widthL = mi - ai;
            widthR = bi - mi;
            if (hmin > 0) && (widthL <= hmin || widthR <= hmin)
                % No se puede subdividir mas de forma segura, aceptar sin correccion
                I = I + (Sleft + Sright);
                n_accept = n_accept + 1;
                accepted(end+1, :) = [ai, mi, bi];
                err_loc(end+1, 1) = err; % se registra el err local (no corregido)
            else
                % Subdividir: push derecha y push izquierda (repartiendo tol/2)
                % Derecha: [mi, di, bi, fmi, fdi, fbi, Sright, tol_i/2]
                stack(end+1, :) = [mi, di, bi, fmi, fdi, fbi, Sright, tol_i / 2];
                % Izquierda: [ai, ci, mi, fai, fci, fmi, Sleft, tol_i/2]
                stack(end+1, :) = [ai, ci, mi, fai, fci, fmi, Sleft, tol_i / 2];
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

