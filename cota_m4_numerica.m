function [M4, info] = cota_M4_numerica(f, a, b, N, safety)
% cota_M4_numerica  Estima una cota superior para |f''''(x)| en [a,b].
% Inputs:
%   f       handle a funcion escalar, ej: @(x) sin(x)
%   a, b    limites con b > a
%   N       (opcional) cantidad de puntos internos para muestrear (default 400)
%   safety  (opcional) factor de seguridad >=1, multiplica el maximo (default 3)
% Outputs:
%   M4      cota superior estimada para |f''''(x)| en [a,b]
%   info    struct con datos de auditoria:
%           .hfd    paso del stencil para derivada cuarta
%           .xs     puntos donde se estimo f''''(x)
%           .f4abs  valores absolutos estimados de f''''(x)
%           .max_no_safety  maximo antes de multiplicar por safety
%           .fevals cantidad de evaluaciones de f
%
% Metodo:
%   Usa el stencil central de 5 puntos para f''''(x):
%     f4(x) ~= [f(x-2h)-4f(x-h)+6f(x)-4f(x+h)+f(x+2h)] / h^4
%   Se define un paso hfd y se evalua en xs dentro de [a+2hfd, b-2hfd].
%
% Notas:
%   - No usa toolboxes. ASCII puro.
%   - Si f es ruidosa o muy curva, aumente N o safety.

    if nargin < 3
        error('Se requieren f, a, b.');
    end
    if ~isa(f, 'function_handle')
        error('f debe ser un function handle.');
    end
    if ~isscalar(a) || ~isscalar(b) || ~(b > a)
        error('Se requiere a < b.');
    end
    if nargin < 4 || isempty(N)
        N = 400;
    end
    if nargin < 5 || isempty(safety)
        safety = 3;
    end
    N = floor(N);
    if N < 10
        N = 10;
    end
    if safety < 1
        safety = 1;
    end

    L = b - a;

    % Paso del stencil: proporcional al tamano de intervalo
    % Elegimos hfd de modo que existan al menos N puntos internos
    % xs en [a+2hfd, b-2hfd]. Tomamos hfd = L / (8*N).
    hfd = L / (8 * N);
    if hfd <= 0
        error('hfd no valido.');
    end

    left  = a + 2 * hfd;
    right = b - 2 * hfd;
    if ~(right > left)
        % Si el intervalo es muy corto, reduzca hfd
        hfd = L / 1000;
        left  = a + 2 * hfd;
        right = b - 2 * hfd;
        if ~(right > left)
            error('Intervalo demasiado corto para stencil de 4ta derivada.');
        end
    end

    % Puntos internos para estimar f''''(x)
    xs = linspace(left, right, N).';

    % Evaluaciones de f para el stencil (vectorizadas por bloques)
    xm2 = xs - 2 * hfd;
    xm1 = xs - 1 * hfd;
    xp1 = xs + 1 * hfd;
    xp2 = xs + 2 * hfd;

    f_m2 = feval(f, xm2);
    f_m1 = feval(f, xm1);
    f_0  = feval(f, xs);
    f_p1 = feval(f, xp1);
    f_p2 = feval(f, xp2);

    if any(~isfinite(f_m2)) || any(~isfinite(f_m1)) || any(~isfinite(f_0)) || any(~isfinite(f_p1)) || any(~isfinite(f_p2))
        error('Valores no finitos de f durante la estimacion de f''''''''.');
    end

    f4_est = (f_m2 - 4 * f_m1 + 6 * f_0 - 4 * f_p1 + f_p2) ./ (hfd ^ 4);
    f4abs = abs(f4_est);

    max_no_safety = max(f4abs);
    M4 = safety * max_no_safety;

    info = struct();
    info.hfd = hfd;
    info.xs = xs;
    info.f4abs = f4abs;
    info.max_no_safety = max_no_safety;
    info.fevals = numel(xs) * 5 + numel(xs);
end

