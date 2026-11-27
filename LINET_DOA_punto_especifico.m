function LINET_Colombia_DOA_punto
% ==========================================================================
% Error absoluto de geolocalizacion por DOA (acimut) para UN SOLO PUNTO
% Solver: Gauss-Newton sobre residuos angulares.
% Muestra el proceso paso a paso de Gauss-Newton.
% ==========================================================================
    %% ================= PARAMETROS =================
    rng(1);
    sigma_theta = deg2rad(1.5);   % incertidumbre angular [rad]
    n_mc        = 3;              % repeticiones Monte Carlo
    mostrar_GN  = true;           % imprime cada iteracion de Gauss-Newton
    detalle_GN  = true;           % imprime tambien residuos, Jacobiano y paso dP

    % ---- PUNTO A EVALUAR (EDITA AQUI) ----
    lat_p = 4.60;                 % latitud  [grados]
    lon_p = -74.08;               % longitud [grados]

    %% ================= SENSORES =================
    lon = [-77.2703; -75.2686; -71.4211];
    lat = [1.2832; 2.9219; 3.8063];
    proj = projcrs(9377);
    [x_s, y_s] = projfwd(proj, lat, lon);
    Sen_pos = [x_s(:), y_s(:)];
    N = size(Sen_pos,1);
fprintf('distancia sensores en metros[grados]: %s\n', x_s(:), y_s(:));
    %% ================= PUNTO REAL =================
    [xp, yp] = projfwd(proj, lat_p, lon_p);
    P_real = [xp, yp];

    %% ================= PROCESAMIENTO =================
    theta_true = atan2(P_real(2) - Sen_pos(:,2), ...
                       P_real(1) - Sen_pos(:,1));

    P_est = zeros(n_mc, 2);
    err   = zeros(n_mc, 1);
    hist_all = cell(n_mc, 1);

    fprintf('\nPunto real: lat = %.5f, lon = %.5f  (x = %.2f, y = %.2f)\n', ...
            lat_p, lon_p, P_real(1), P_real(2));
    fprintf('Sigma angular: %.3f grados | Repeticiones: %d\n', rad2deg(sigma_theta), n_mc);
    fprintf('Acimuts verdaderos [grados]: %s\n', mat2str(rad2deg(theta_true)', 5));

    for m = 1:n_mc
        theta_meas = theta_true + randn(N,1)*sigma_theta;

        fprintf('\n==========================================================\n');
        fprintf(' REPETICION %d\n', m);
        fprintf('==========================================================\n');
        fprintf('Acimuts medidos  [grados]: %s\n', mat2str(rad2deg(theta_meas)', 5));
        fprintf('Ruido agregado   [grados]: %s\n', mat2str(rad2deg(theta_meas - theta_true)', 4));

        [P_doa, n_it, hist] = solve_DOA_GN(Sen_pos, theta_meas, mostrar_GN, detalle_GN, P_real);

        P_est(m,:)   = P_doa;
        err(m)       = norm(P_real - P_doa);
        hist_all{m}  = hist;
        [lat_e, lon_e] = projinv(proj, P_doa(1), P_doa(2));

        fprintf('\n>> Resultado: lat = %.5f, lon = %.5f | Error = %.1f m | Iteraciones GN = %d\n', ...
                lat_e, lon_e, err(m), n_it);
    end

    %% ================= RESUMEN =================
    fprintf('\n----------------------------------------------------------\n');
    fprintf('Resumen Monte Carlo (%d repeticiones):\n', n_mc);
    fprintf('  Error medio    : %.1f m\n', mean(err));
    fprintf('  Mediana        : %.1f m\n', median(err));
    fprintf('  Desv. estandar : %.1f m\n', std(err));
    fprintf('  Error minimo   : %.1f m\n', min(err));
    fprintf('  Error maximo   : %.1f m\n', max(err));

    %% ================= FIGURA 1: ESTIMACIONES =================
    figure('Color','w','Name','DOA - Punto unico','Position',[60 80 700 550]);
    hold on
    plot(P_est(:,1), P_est(:,2), '.', 'Color',[0.4 0.4 0.4], 'MarkerSize',14);
    plot(x_s, y_s, '^', 'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'MarkerSize',8);
    plot(P_real(1), P_real(2), 'rp', 'MarkerFaceColor','r', 'MarkerSize',12);
    plot(mean(P_est(:,1)), mean(P_est(:,2)), 'bx', 'MarkerSize',10, 'LineWidth',2);
    axis equal; grid on; box on
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('DOA - Error medio: %.1f m', mean(err)));
    legend('Estimaciones','Sensores','Punto real','Media estimada','Location','best');
    hold off

    %% ================= FIGURA 2: TRAYECTORIA DE GAUSS-NEWTON =================
    figure('Color','w','Name','DOA - Convergencia Gauss-Newton','Position',[780 80 700 550]);
    hold on
    cols = lines(n_mc);
    for m = 1:n_mc
        h = hist_all{m};
        plot(h(:,1), h(:,2), '-o', 'Color',cols(m,:), 'MarkerSize',5, ...
             'MarkerFaceColor',cols(m,:), 'DisplayName', sprintf('Repeticion %d', m));
        plot(h(1,1), h(1,2), 's', 'Color',cols(m,:), 'MarkerSize',9, ...
             'HandleVisibility','off');            % cuadrado = semilla
    end
    plot(P_real(1), P_real(2), 'rp', 'MarkerFaceColor','r', 'MarkerSize',14, ...
         'DisplayName','Punto real');
    axis equal; grid on; box on
    xlabel('x (m)'); ylabel('y (m)');
    title('Trayectoria de Gauss-Newton (cuadrado = semilla)');
    legend('Location','best');
    hold off
end

%% ==========================================================================
%  DOA - Gauss-Newton (incognitas x,y) sobre residuos angulares
%  Modelo:    theta_i(P) = atan2(Py - Sy_i, Px - Sx_i)
%  Actualiz.: dP = (J'J)^-1 J'r ,  P <- P + dP
%  Salidas: P = estimacion final, n_it = iteraciones, hist = trayectoria [x y]
%% ==========================================================================
function [P, n_it, hist] = solve_DOA_GN(S, theta_meas, verbose, detalle, P_real)
    % Semilla: interseccion lineal de rectas de acimut (minimos cuadrados)
    A = [sin(theta_meas), -cos(theta_meas)];
    b = sin(theta_meas).*S(:,1) - cos(theta_meas).*S(:,2);
    P = (A \ b)';
    hist = P;

    if verbose
        fprintf('\nSemilla (interseccion lineal de rectas):\n');
        fprintf('  x = %.2f m, y = %.2f m | error = %.2f m\n', P(1), P(2), norm(P - P_real));
    end

    n_it = 0;
    for it = 1:20
        dx = P(1) - S(:,1);
        dy = P(2) - S(:,2);
        d2 = dx.^2 + dy.^2;
        r = atan2(sin(theta_meas - atan2(dy,dx)), ...
                  cos(theta_meas - atan2(dy,dx)));   % residuo envuelto a [-pi,pi]
        J = [-dy./d2, dx./d2];
        dP = ((J'*J) + 1e-18*eye(2)) \ (J'*r);

        limitado = false;
        if norm(dP) > 2e5
            dP = dP*(2e5/norm(dP));                  % limite de paso
            limitado = true;
        end
        P = P + dP';
        n_it = it;
        hist(end+1,:) = P; %#ok<AGROW>

        if verbose
            fprintf('\n--- Iteracion %d ---\n', it);
            if detalle
                fprintf('  Residuos r [grados]    : %s\n', mat2str(rad2deg(r)', 4));
                fprintf('  Jacobiano J [rad/m] (una fila por sensor: [dtheta/dx  dtheta/dy]):\n');
                for k = 1:size(J,1)
                    fprintf('     S%d: [% .4e  % .4e]\n', k, J(k,1), J(k,2));
                end
                fprintf('  Paso dP = (J''J)^-1 J''r : [%.3f  %.3f] m', dP(1), dP(2));
                if limitado, fprintf('  (limitado a 2e5 m)'); end
                fprintf('\n');
            end
            fprintf('  |r| = %.3e rad | |dP| = %.4f m\n', norm(r), norm(dP));
            fprintf('  Nueva posicion: x = %.2f m, y = %.2f m | error = %.2f m\n', ...
                    P(1), P(2), norm(P - P_real));
        end

        if norm(dP) < 1e-3
            if verbose, fprintf('\n  Convergio: |dP| < 1e-3 m\n'); end
            break;
        end
    end
end