function LINET_Colombia_DOA_general
% ==========================================================================
% Mapa de calor GENERAL del error absoluto de geolocalizacion por DOA
% (Direction Of Arrival, acimut). Solver: Gauss-Newton sobre residuos
% angulares. Sin croquis de Colombia.
% ==========================================================================

    %% ================= PARAMETROS =================
    rng(1);
    sigma_theta = deg2rad(0.5);   % incertidumbre angular [rad] (~1.5 deg)
    grid_res = 120;               % puntos por eje
    n_mc     = 5;                 % repeticiones Monte Carlo por punto
    margen   = 10000e3;            % margen alrededor de la red [m]
    cmax     = 500000;             % tope de la barra de color [m]

    %% ================= SENSORES =================
    lon = [-77.2703; -75.2686; -71.4211];
    lat = [1.2832; 2.9219; 3.8063];
    proj = projcrs(9377);
    [x_s, y_s] = projfwd(proj, lat, lon);
    Sen_pos = [x_s(:), y_s(:)];
    N = size(Sen_pos,1);

    %% ================= MALLA RECTANGULAR =================
    x_vec = linspace(min(x_s)-margen, max(x_s)+margen, grid_res);
    y_vec = linspace(min(y_s)-margen, max(y_s)+margen, grid_res);
    [X, Y] = meshgrid(x_vec, y_vec);

    Err_DOA = nan(size(X));

    %% ================= PROCESAMIENTO =================
    fprintf('Calculando %d puntos (x%d MC)... ', numel(X), n_mc);
    t_ini = tic;

    for i = 1:numel(X)
        P_real = [X(i), Y(i)];
        theta_true = atan2(P_real(2) - Sen_pos(:,2), ...
                           P_real(1) - Sen_pos(:,1));

        e = zeros(n_mc,1);
        for m = 1:n_mc
            theta_meas = theta_true + randn(N,1)*sigma_theta;
            P_doa = solve_DOA_GN(Sen_pos, theta_meas);
            e(m)  = norm(P_real - P_doa);
        end
        Err_DOA(i) = mean(e);      % error absoluto medio [m]
    end
    fprintf('Hecho (%.1f s).\n', toc(t_ini));

    %% ================= VISUALIZACION =================
    figure('Color','w','Name','Mapa de calor DOA - Red LINET', ...
           'Position',[60 80 800 600]);
    plot_error_map(x_vec, y_vec, Err_DOA, x_s, y_s, ...
        {'DOA','Angulo de llegada (acimut)'}, cmax);

    fprintf('\nError absoluto medio sobre toda la malla:\n');
    fprintf('  DOA : %.1f m\n', mean(Err_DOA(:),'omitnan'));
end

%% ==========================================================================
%  GRAFICADOR: mapa de calor + sensores
%% ==========================================================================
function plot_error_map(x_vec, y_vec, Err, x_s, y_s, titulo, cmax)
    hold on
    imagesc(x_vec, y_vec, Err);
    set(gca, 'YDir','normal');
    colormap(gca, turbo);
    clim([0 cmax]);                    % MATLAB antiguo: caxis([0 cmax])
    cb = colorbar; cb.Label.String = 'Error absoluto [m]';

    plot(x_s, y_s, '^', 'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'MarkerSize',6);

    axis equal tight
    grid on; box on
    xlabel('x (m)'); ylabel('y (m)');
    title(titulo);
    hold off
end

%% ==========================================================================
%  DOA - Gauss-Newton (incognitas x,y) sobre residuos angulares
%  Modelo:   theta_i(P) = atan2(Py - Sy_i, Px - Sx_i)
%  Jacobiano: d theta_i / dP = [ -dy_i/d_i^2 , dx_i/d_i^2 ]
%% ==========================================================================
function P = solve_DOA_GN(S, theta_meas)
    % Semilla: interseccion lineal de rectas de acimut (minimos cuadrados)
    A = [sin(theta_meas), -cos(theta_meas)];
    b = sin(theta_meas).*S(:,1) - cos(theta_meas).*S(:,2);
    P = (A \ b)';

    for it = 1:20
        dx = P(1) - S(:,1);
        dy = P(2) - S(:,2);
        d2 = dx.^2 + dy.^2;

        r = atan2(sin(theta_meas - atan2(dy,dx)), ...
                  cos(theta_meas - atan2(dy,dx)));   % residuo envuelto a [-pi,pi]
        J = [-dy./d2, dx./d2];

        dP = ((J'*J) + 1e-18*eye(2)) \ (J'*r);
        if norm(dP) > 2e5, dP = dP*(2e5/norm(dP)); end   % limite de paso
        P = P + dP';
        if norm(dP) < 1e-3, break; end
    end
end