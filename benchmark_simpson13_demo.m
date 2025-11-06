function results = benchmark_simpson13_demo(tol)
% benchmark_simpson13_demo  Compara Simpson 1/3 auto vs adaptativo en 2 casos.
% Input:
%   tol  tolerancia absoluta objetivo (default 1e-8)
% Output:
%   results  struct con resultados para sin(x) y funcion con pico
%
% Requiere en el path:
%   - simpson13_f_auto.m
%   - simpson13_adaptativo_f.m
%
% Uso rapido:
%   results = benchmark_simpson13_demo(1e-8);

    if nargin < 1
        tol = 1e-8;
    end

    % ==== Caso A: f(x) = sin(x) en [0, pi]  ====
    fA = @(x) sin(x);
    aA = 0;
    bA = pi;
    IexA = 2.0;

    [IA_auto, hA, nA, errA_est, feA_auto] = simpson13_f_auto(fA, aA, bA, tol);
    [IA_adap, infoA] = simpson13_adaptativo_f(fA, aA, bA, tol);

    errA_abs_auto = abs(IA_auto - IexA);
    errA_abs_adap = abs(IA_adap - IexA);

    % ==== Caso B: g(x) = 1/(1+100(x-0.7)^2) en [0,1] ====
    gB = @(x) 1.0 ./ (1.0 + 100.0 * (x - 0.7) .* (x - 0.7));
    aB = 0.0;
    bB = 1.0;
    % Integral exacta: (1/10) * [atan(10(x-0.7))]_{0}^{1} = (atan(3)+atan(7))/10
    IexB = (atan(3.0) + atan(7.0)) / 10.0;

    [IB_auto, hB, nB, errB_est, feB_auto] = simpson13_f_auto(gB, aB, bB, tol);
    [IB_adap, infoB] = simpson13_adaptativo_f(gB, aB, bB, tol);

    errB_abs_auto = abs(IB_auto - IexB);
    errB_abs_adap = abs(IB_adap - IexB);

    % ==== Imprimir resumen ====
    fprintf('================= Benchmark Simpson 1/3 (tol = %.3e) =================\n', tol);
    fprintf('Caso A: f(x) = sin(x) en [0, pi], I_exacta = 2\n');
    fprintf(' Metodo        I_aprox            Error_abs        Eval_f     Extra\n');
    fprintf(' Auto      % .12e   % .3e   %8d     n=%d, h=%.3e, est=%.3e\n', ...
        IA_auto, errA_abs_auto, feA_auto, nA, hA, errA_est);
    fprintf(' Adapt     % .12e   % .3e   %8d     acept=%d, split=%d\n', ...
        IA_adap, errA_abs_adap, infoA.fevals, infoA.n_accept, infoA.n_split);

    fprintf('\nCaso B: g(x) = 1/(1+100(x-0.7)^2) en [0,1], I_exacta = (atan(3)+atan(7))/10\n');
    fprintf(' Metodo        I_aprox            Error_abs        Eval_f     Extra\n');
    fprintf(' Auto      % .12e   % .3e   %8d     n=%d, h=%.3e, est=%.3e\n', ...
        IB_auto, errB_abs_auto, feB_auto, nB, hB, errB_est);
    fprintf(' Adapt     % .12e   % .3e   %8d     acept=%d, split=%d\n', ...
        IB_adap, errB_abs_adap, infoB.fevals, infoB.n_accept, infoB.n_split);
    fprintf('======================================================================\n');

    % ==== Devolver struct con resultados ====
    A = struct();
    A.exacto = IexA;
    A.auto = struct('I', IA_auto, 'h', hA, 'n', nA, 'err_abs', errA_abs_auto, 'err_est', errA_est, 'fevals', feA_auto);
    A.adapt = struct('I', IA_adap, 'err_abs', errA_abs_adap, 'fevals', infoA.fevals, 'n_accept', infoA.n_accept, 'n_split', infoA.n_split);

    B = struct();
    B.exacto = IexB;
    B.auto = struct('I', IB_auto, 'h', hB, 'n', nB, 'err_abs', errB_abs_auto, 'err_est', errB_est, 'fevals', feB_auto);
    B.adapt = struct('I', IB_adap, 'err_abs', errB_abs_adap, 'fevals', infoB.fevals, 'n_accept', infoB.n_accept, 'n_split', infoB.n_split);

    results = struct();
    results.casoA = A;
    results.casoB = B;
end

