function LINET_Colombia_TOA_TDOA_punto
% ==========================================================================
% Error absoluto de geolocalizacion (TDOA y TOA) para UN SOLO PUNTO
% Red de sensores tipo LINET. Solver: Gauss-Newton.
% Monte Carlo + impresion de cada iteracion de Gauss-Newton.
% ==========================================================================

    %% ================= PARAMETROS =================
    v_luz    = 299792458;     % [m/s]
    rng(1);
    sigma_t  = 50e-9;         % incertidumbre temporal [s]
    n_mc     = 3;             % repeticiones Monte Carlo
    mostrar_GN = true;        % imprime cada iteracion de Gauss-Newton
    detalle_GN = true;        % imprime tambien residuos, Jacobiano y paso

    % ---- PUNTO A EVALUAR (EDITA AQUI) ----
    lat_p = 4.60;             % latitud  [grados]
    lon_p = -74.08;           % longitud [grados]

    %% ================= SENSORES =================
    lon = [-77.2703; -75.2686; -71.4211];
    lat = [1.2832; 2.9219; 3.8063];
    proj = projcrs(9377);
    [x_s, y_s] = projfwd(proj, lat, lon);
    Sen_pos = [x_s(:), y_s(:)];
    N = size(Sen_pos,1);

    %% ================= PUNTO REAL =================
    [xp, yp] = projfwd(proj, lat_p, lon_p);
    P_real = [xp, yp];
    d = sqrt(sum((Sen_pos - P_real).^2, 2));
    P0 = mean(Sen_pos, 1);                 % semilla comun punto donde de va a iniciar

    fprintf('\nPunto real: lat = %.5f, lon = %.5f  (x = %.2f, y = %.2f)\n', ...
            lat_p, lon_p, P_real(1), P_real(2));
    fprintf('Sigma temporal: %.1f ns | Repeticiones: %d\n', sigma_t*1e9, n_mc);
    fprintf('Semilla P0 (centroide de sensores): x = %.2f, y = %.2f | error = %.1f m\n', ...
            P0(1), P0(2), norm(P0 - P_real));
    fprintf('Distancias reales a sensores [km]: %s\n', mat2str(d'/1e3, 6));
    fprintf('Tiempos de llegada verdaderos [ms]: %s\n', mat2str(d'/v_luz*1e3, 8));

    P_tdoa_all = zeros(n_mc,2);  P_toa_all = zeros(n_mc,2);
    e_tdoa = zeros(n_mc,1);      e_toa = zeros(n_mc,1);
    h_tdoa = cell(n_mc,1);       h_toa = cell(n_mc,1);

    %% ================= PROCESAMIENTO =================
    for m = 1:n_mc
        ruido    = randn(N,1)*sigma_t;
        TOA_meas = d./v_luz + ruido;

        fprintf('\n##########################################################\n');
        fprintf(' REPETICION %d\n', m);
        fprintf('##########################################################\n');
        fprintf('Ruido agregado [ns]: %s\n', mat2str(ruido'*1e9, 4));

        % ---------- TDOA ----------
        fprintf('\n=================== TDOA (ref. sensor 1) ===================\n');
        [P_tdoa, n1, h1] = solve_TDOA_GN(Sen_pos, TOA_meas, v_luz, P0, ...
                                         mostrar_GN, detalle_GN, P_real);
        % ---------- TOA ----------
        fprintf('\n=================== TOA (t0 desconocido) ===================\n');
        [P_toa, n2, h2] = solve_TOA_GN(Sen_pos, TOA_meas, v_luz, P0, ...
                                       mostrar_GN, detalle_GN, P_real);

        P_tdoa_all(m,:) = P_tdoa;      P_toa_all(m,:) = P_toa;
        e_tdoa(m) = norm(P_real - P_tdoa);
        e_toa(m)  = norm(P_real - P_toa);
        h_tdoa{m} = h1;                h_toa{m} = h2;

        [la1, lo1] = projinv(proj, P_tdoa(1), P_tdoa(2));
        [la2, lo2] = projinv(proj, P_toa(1),  P_toa(2));
        fprintf('\n>> TDOA: lat = %.5f, lon = %.5f | Error = %.1f m | Iter GN = %d\n', la1, lo1, e_tdoa(m), n1);
        fprintf('>> TOA : lat = %.5f, lon = %.5f | Error = %.1f m | Iter GN = %d\n',  la2, lo2, e_toa(m),  n2);
    end

    %% ================= RESUMEN =================
    fprintf('\n----------------------------------------------------------\n');
    fprintf('Resumen Monte Carlo (%d repeticiones):\n', n_mc);
    fprintf('            %12s %12s\n', 'TDOA', 'TOA');
    fprintf('  Medio   : %12.1f %12.1f m\n', mean(e_tdoa),   mean(e_toa));
    fprintf('  Mediana : %12.1f %12.1f m\n', median(e_tdoa), median(e_toa));
    fprintf('  Desv.est: %12.1f %12.1f m\n', std(e_tdoa),    std(e_toa));
    fprintf('  Minimo  : %12.1f %12.1f m\n', min(e_tdoa),    min(e_toa));
    fprintf('  Maximo  : %12.1f %12.1f m\n', max(e_tdoa),    max(e_toa));

    %% ================= FIGURA 1: ESTIMACIONES =================
    figure('Color','w','Name','TDOA / TOA - Estimaciones','Position',[60 80 1300 520]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    nexttile
    plot_estimaciones(Sen_pos, P_real, P_tdoa_all, e_tdoa, 'TDOA');
    nexttile
    plot_estimaciones(Sen_pos, P_real, P_toa_all,  e_toa,  'TOA');

    %% ================= FIGURA 2: TRAYECTORIA DE GAUSS-NEWTON =================
    figure('Color','w','Name','TDOA / TOA - Convergencia Gauss-Newton','Position',[60 80 1300 520]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
    nexttile
    plot_trayectoria(Sen_pos, P_real, P0, h_tdoa, 'TDOA');
    nexttile
    plot_trayectoria(Sen_pos, P_real, P0, h_toa,  'TOA');
end

%% ==========================================================================
%  GRAFICADORES
%% ==========================================================================
function plot_estimaciones(S, P_real, P_est, err, nombre)
    hold on
    plot(P_est(:,1), P_est(:,2), '.', 'Color',[0.4 0.4 0.4], 'MarkerSize',14);
    plot(S(:,1), S(:,2), '^', 'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'MarkerSize',8);
    plot(P_real(1), P_real(2), 'rp', 'MarkerFaceColor','r', 'MarkerSize',12);
    plot(mean(P_est(:,1)), mean(P_est(:,2)), 'bx', 'MarkerSize',10, 'LineWidth',2);
    axis equal; grid on; box on
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('%s - Error medio: %.1f m', nombre, mean(err)));
    legend('Estimaciones','Sensores','Punto real','Media estimada','Location','best');
    hold off
end

function plot_trayectoria(S, P_real, P0, hist_all, nombre)
    hold on
    n_mc = numel(hist_all);
    cols = lines(n_mc);
    for m = 1:n_mc
        h = hist_all{m};
        plot(h(:,1), h(:,2), '-o', 'Color',cols(m,:), 'MarkerSize',5, ...
             'MarkerFaceColor',cols(m,:), 'DisplayName', sprintf('Repeticion %d', m));
    end
    plot(P0(1), P0(2), 'ks', 'MarkerSize',10, 'DisplayName','Semilla P0');
    plot(S(:,1), S(:,2), '^', 'MarkerFaceColor','w', 'MarkerEdgeColor','k', ...
         'MarkerSize',8, 'DisplayName','Sensores');
    plot(P_real(1), P_real(2), 'rp', 'MarkerFaceColor','r', 'MarkerSize',14, ...
         'DisplayName','Punto real');
    axis equal; grid on; box on
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('%s - Trayectoria de Gauss-Newton', nombre));
    legend('Location','best');
    hold off
end

%% ==========================================================================
%  TDOA - Gauss-Newton (incognitas x,y). Referencia: sensor 1
%  Residuo:  r_i = c*(t_i - t_1) - (d_i - d_1),  i = 2..N
%  Salidas: P, n_it, hist (trayectoria [x y])
%% ==========================================================================
function [P, n_it, hist] = solve_TDOA_GN(S, TOA, c, P0, verbose, detalle, P_real)
    N  = size(S,1);
    dt = TOA - TOA(1);
    idx = 2:N;
    P = P0(:)';
    hist = P;

    if verbose
        fprintf('Diferencias de tiempo medidas dt_i1 [ns]: %s\n', mat2str(dt(idx)'*1e9, 6));
        fprintf('Inicio: x = %.2f m, y = %.2f m | error = %.2f m\n', P(1), P(2), norm(P - P_real));
    end

    n_it = 0;
    for it = 1:20
        dist = sqrt(sum((P - S).^2, 2));
        r  = c*dt(idx) - (dist(idx) - dist(1));
        u  = (P - S) ./ dist;
        Jf = u(idx,:) - u(1,:);
        dP = ((Jf'*Jf) + 1e-9*eye(2)) \ (Jf'*r);

        limitado = false;
        if norm(dP) > 2e5
            dP = dP*(2e5/norm(dP));
            limitado = true;
        end
        P = P + dP';
        n_it = it;
        hist(end+1,:) = P; %#ok<AGROW>

        if verbose
            fprintf('\n--- TDOA Iteracion %d ---\n', it);
            if detalle
                fprintf('  Distancias actuales [km]: %s\n', mat2str(dist'/1e3, 6));
                fprintf('  Residuos r [m]          : %s\n', mat2str(r', 6));
                fprintf('  Jacobiano J (fila = par i-1: [dr/dx  dr/dy]):\n');
                for k = 1:size(Jf,1)
                    fprintf('     S%d-S1: [% .5f  % .5f]\n', idx(k), Jf(k,1), Jf(k,2));
                end
                fprintf('  Paso dP = (J''J + 1e-9 I)^-1 J''r : [%.3f  %.3f] m', dP(1), dP(2));
                if limitado, fprintf('  (limitado a 2e5 m)'); end
                fprintf('\n');
            end
            fprintf('  |r| = %.3e m | |dP| = %.4f m\n', norm(r), norm(dP));
            fprintf('  Nueva posicion: x = %.2f m, y = %.2f m | error = %.2f m\n', ...
                    P(1), P(2), norm(P - P_real));
        end

        if norm(dP) < 1e-3
            if verbose, fprintf('\n  Convergio: |dP| < 1e-3 m\n'); end
            break;
        end
    end
end

%% ==========================================================================
%  TOA - Gauss-Newton con incognitas [x, y, b], donde b = c*t0
%  Residuo:  r_i = c*t_i - (d_i + b)
%  Salidas: P, n_it, hist (trayectoria [x y])
%% ==========================================================================
function [P, n_it, hist] = solve_TOA_GN(S, TOA, c, P0, verbose, detalle, P_real)
    N = size(S,1);
    z = c*TOA;
    d0 = sqrt(sum((P0 - S).^2, 2));
    p = [P0(:); mean(z - d0)];
    hist = p(1:2)';

    if verbose
        fprintf('Distancias medidas c*t_i [km]: %s\n', mat2str(z'/1e3, 8));
        fprintf('Inicio: x = %.2f m, y = %.2f m, b = c*t0 = %.2f m | error = %.2f m\n', ...
                p(1), p(2), p(3), norm(p(1:2)' - P_real));
    end

    n_it = 0;
    for it = 1:20
        dvec = p(1:2)' - S;
        dist = sqrt(sum(dvec.^2, 2));
        r = z - (dist + p(3));
        J = [dvec./dist, ones(N,1)];
        dp = (J'*J + 1e-9*eye(3)) \ (J'*r);

        limitado = false;
        if norm(dp(1:2)) > 2e5
            dp = dp*(2e5/norm(dp(1:2)));
            limitado = true;
        end
        p = p + dp;
        n_it = it;
        hist(end+1,:) = p(1:2)'; %#ok<AGROW>

        if verbose
            fprintf('\n--- TOA Iteracion %d ---\n', it);
            if detalle
                fprintf('  Distancias actuales [km]: %s\n', mat2str(dist'/1e3, 6));
                fprintf('  Residuos r [m]          : %s\n', mat2str(r', 6));
                fprintf('  Jacobiano J (fila = sensor: [dr/dx  dr/dy  dr/db]):\n');
                for k = 1:N
                    fprintf('     S%d: [% .5f  % .5f  % .5f]\n', k, J(k,1), J(k,2), J(k,3));
                end
                fprintf('  Paso dp : [dx = %.3f, dy = %.3f, db = %.3f] m', dp(1), dp(2), dp(3));
                if limitado, fprintf('  (limitado a 2e5 m)'); end
                fprintf('\n');
            end
            fprintf('  |r| = %.3e m | |dp(x,y)| = %.4f m\n', norm(r), norm(dp(1:2)));
            fprintf('  Nueva posicion: x = %.2f m, y = %.2f m, b = %.2f m | error = %.2f m\n', ...
                    p(1), p(2), p(3), norm(p(1:2)' - P_real));
        end

        if norm(dp(1:2)) < 1e-3
            if verbose, fprintf('\n  Convergio: |dp(x,y)| < 1e-3 m\n'); end
            break;
        end
    end
    P = p(1:2)';
end