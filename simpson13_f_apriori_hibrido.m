function [I, h, n, err_est, M4] = simpson13_f_apriori_hibrido(f, a, b, tol, N_M4, safety, kmax)
% simpson13_f_apriori_hibrido
%   Simpson 1/3 compuesto con seleccion a priori de n usando una cota
%   numerica de |f''''(x)|, y verificacion a posteriori por mallas anidadas.
%
% Inputs:
%   f       handle a funcion escalar
%   a, b    limites con b > a
%   tol     tolerancia absoluta deseada
%   N_M4    (opcional) puntos para estimar M4 (default 400)
%   safety  (opcional) factor de seguridad para M4 (default 3)
%   kmax    (opcional) maximo numero de refinamientos a posteriori (default 10)
%
% Outputs:
%   I       aproximacion final de la integral
%   h       paso final de malla (h = (b - a)/n)
%   n       numero final de subintervalos (par)
%   err_est estimacion a posteriori |I_{2n} - I_n| / 15
%   M4      cota usada para |f''''(x)|
%
% Requisitos:
%   - Requiere tener en el path el archivo cota_M4_numerica.m
%   - No usa toolboxes. ASCII puro.
%
% Metodo:
%   1) Estima M4 numericamente (central 5 puntos) con factor de seguridad.
%   2) Fija hreq = ((180*tol)/((b-a)*M4))^(1/4) y n = ceil((b-a)/hreq),
%      ajustando n al par inmediato y n>=2.
%   3) Calcula I_n por Simpson 1/3. Luego I_2n y err_est = |I_2n - I_n|/15.
%   4) Si err_est > tol, duplica n hasta cumplir o agotar kmax.
%
% Nota:
%   - Si M4 resulta 0 (p.ej. polinomio grado <= 3), usa n = 2 y valida.

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
    if nargin < 5 || isempty(N_M4)
        N_M4 = 400;
    end
    if nargin < 6 || isempty(safety)
        safety = 3;
    end
    if nargin < 7 || isempty(kmax)
        kmax = 10;
    end

    L = b - a;

    % 1) Estimar M4 numericamente
    [M4, ~] = cota_M4_numerica(f, a, b, N_M4, safety);

    % 2) Elegir n a priori con la cota
    if M4 == 0
        n = 2;
    else
        hreq = (180 * tol / (L * M4))^(1/4);
        if ~(isfinite(hreq)) || hreq <= 0
            error('No se pudo calcular h requerido. Revise M4 y tol.');
        end
        n = ceil(L / hreq);
        if n < 2
            n = 2;
        end
    end
    if mod(n, 2) ~= 0
        n = n + 1;
    end

    % 3) Simpson 1/3 con verificacion a posteriori
    I_n = simpson_once_local(f, a, b, n);
    n2 = 2 * n;
    I_n2 = simpson_once_local(f, a, b, n2);
    err_est = abs(I_n2 - I_n) / 15;

    kk = 0;
    while err_est > tol && kk < kmax
        n = n2;
        I_n = I_n2;
        n2 = 2 * n;
        I_n2 = simpson_once_local(f, a, b, n2);
        err_est = abs(I_n2 - I_n) / 15;
        kk = kk + 1;
    end

    I = I_n2;
    n = n2;
    h = L / n;
end

function S = simpson_once_local(f, a, b, n)
% simpson_once_local  Simpson 1/3 compuesto en malla uniforme con n par.

    if ~isscalar(n) || n < 2 || mod(n, 2) ~= 0 || floor(n) ~= n
        error('Simpson requiere n par, n >= 2.');
    end

    h = (b - a) / n;

    fa = feval(f, a);
    fb = feval(f, b);
    if ~isfinite(fa) || ~isfinite(fb)
        error('Valores no finitos de f en los extremos.');
    end

    s4 = 0;
    i = 1;
    while i <= n - 1
        xi = a + i * h;
        fi = feval(f, xi);
        if ~isfinite(fi)
            error('Valor no finito de f en un nodo (impar).');
        end
        s4 = s4 + fi;
        i = i + 2;
    end

    s2 = 0;
    i = 2;
    while i <= n - 2
        xi = a + i * h;
        fi = feval(f, xi);
        if ~isfinite(fi)
            error('Valor no finito de f en un nodo (par).');
        end
        s2 = s2 + fi;
        i = i + 2;
    end

    S = (h / 3) * (fa + fb + 4 * s4 + 2 * s2);
end

