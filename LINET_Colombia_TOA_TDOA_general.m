function LINET_Colombia_TOA_TDOA_general
% ==========================================================================
% Mapa de calor GENERAL del error absoluto de geolocalizacion (TDOA y TOA)
% Red simulada de 17 sensores tipo LINET. Sin croquis de Colombia:
% la malla cubre todo el rectangulo que rodea la red + un margen.
% ==========================================================================

    %% ================= PARAMETROS =================
    v_luz    = 299792458;     % [m/s]
    rng(1);
    sigma_t  = 50e-9;         % incertidumbre temporal [s]
    grid_res = 120;           % puntos por eje
    n_mc     = 5;            % repeticiones Monte Carlo por punto
    margen   = 10000e3;         % margen alrededor de la red [m]
    cmax     = 50000;          % tope de la barra de color [m]

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

    Err_TDOA = nan(size(X));
    Err_TOA  = nan(size(X));

    %% ================= PROCESAMIENTO =================
    fprintf('Calculando %d puntos (x%d MC)... ', numel(X), n_mc);
    t_ini = tic;
    P0 = mean(Sen_pos, 1);

    for i = 1:numel(X)
        P_real = [X(i), Y(i)];
        d = sqrt(sum((Sen_pos - P_real).^2, 2));

        e1 = zeros(n_mc,1); e2 = e1;
        for m = 1:n_mc
            TOA_meas = d./v_luz + randn(N,1)*sigma_t;

            P_tdoa = solve_TDOA_GN(Sen_pos, TOA_meas, v_luz, P0);
            P_toa  = solve_TOA_GN(Sen_pos, TOA_meas, v_luz, P0);

            e1(m) = norm(P_real - P_tdoa);
            e2(m) = norm(P_real - P_toa);
        end
        Err_TDOA(i) = mean(e1);    % error absoluto medio [m]
        Err_TOA(i)  = mean(e2);
    end
    fprintf('Hecho (%.1f s).\n', toc(t_ini));

    %% ================= VISUALIZACION =================
    figure('Color','w','Name','Mapas de calor generales - Red LINET', ...
           'Position',[60 80 1300 560]);
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile
    plot_error_map(x_vec, y_vec, Err_TDOA, x_s, y_s, ...
        {'TDOA','Diferencia de tiempos de llegada'}, cmax);
    nexttile
    plot_error_map(x_vec, y_vec, Err_TOA, x_s, y_s, ...
        {'TOA','Tiempo absoluto (t0 desconocido)'}, cmax);

    sgtitle('Error absoluto de geolocalizacion de rayos - 17 sensores', ...
            'FontSize',14,'FontWeight','bold');

    fprintf('\nError absoluto medio sobre toda la malla:\n');
    fprintf('  TDOA : %.1f m\n', mean(Err_TDOA(:),'omitnan'));
    fprintf('  TOA  : %.1f m\n', mean(Err_TOA(:), 'omitnan'));
end

%% ==========================================================================
%  GRAFICADOR: mapa de calor + sensores (sin contorno)
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
%  TDOA - Gauss-Newton (incognitas x,y). Referencia: sensor 1
%% ==========================================================================
function P = solve_TDOA_GN(S, TOA, c, P0)
    N  = size(S,1);
    dt = TOA - TOA(1);
    idx = 2:N;
    P = P0(:)';
    for it = 1:20
        dist = sqrt(sum((P - S).^2, 2));
        r  = c*dt(idx) - (dist(idx) - dist(1));
        u  = (P - S) ./ dist;
        Jf = u(idx,:) - u(1,:);
        dP = ((Jf'*Jf) + 1e-9*eye(2)) \ (Jf'*r);
        if norm(dP) > 2e5, dP = dP*(2e5/norm(dP)); end
        P = P + dP';
        if norm(dP) < 1e-3, break; end
    end
end

%% ==========================================================================
%  TOA - Gauss-Newton con incognitas [x, y, b], donde b = c*t0
%% ==========================================================================
function P = solve_TOA_GN(S, TOA, c, P0)
    N = size(S,1);
    z = c*TOA;
    d0 = sqrt(sum((P0 - S).^2, 2));
    p = [P0(:); mean(z - d0)];
    for it = 1:20
        dvec = p(1:2)' - S;
        dist = sqrt(sum(dvec.^2, 2));
        r = z - (dist + p(3));
        J = [dvec./dist, ones(N,1)];
        dp = (J'*J + 1e-9*eye(3)) \ (J'*r);
        if norm(dp(1:2)) > 2e5, dp = dp*(2e5/norm(dp(1:2))); end
        p = p + dp;
        if norm(dp(1:2)) < 1e-3, break; end
    end
    P = p(1:2)';
end