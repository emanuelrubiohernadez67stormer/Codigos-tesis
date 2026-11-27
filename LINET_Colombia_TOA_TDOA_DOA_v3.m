function LINET_Colombia_TOA_TDOA_DOA_v3
% ==========================================================================
% Mapa de calor del error de geolocalización (TDOA, TOA, DOA) sobre Colombia
% Red simulada de 17 sensores tipo LINET. Mapas en proyección MAGNA-SIRGAS
% Origen Nacional (EPSG:9377), con el croquis de Colombia debajo del mapa.
% ==========================================================================

    %% ================= PARÁMETROS GENERALES =================
    v_luz       = 299792458;          % [m/s]
    rng(1);
    sigma_t     = 100e-9;              % incertidumbre temporal [s]                                  
    grid_res    = 120;                % puntos por eje (subir = más detalle)
    n_mc        = 5;                  % repeticiones Monte Carlo por punto

    %% ================= SENSORES =================
    lon = [-77.2703; -75.2686; -71.4211; -76.2824; -74.0727; -75.2166; -73.3448; -74.6447; -73.8388; -75.5545; -70.7452; -72.4610; -73.3188; -72.6169; -74.8266; -75.8665; -76.6724];
    lat = [1.2832; 2.9219; 3.8063; 3.5202; 4.6387; 4.4306; 5.5751; 5.8612; 7.0578; 6.2774; 7.0317; 7.8381; 9.5549; 11.1676; 10.9855; 8.7745; 7.6820];

    %% ================= PROYECCIÓN Y CONTORNO =================
    proj = projcrs(9377);
    [x_s, y_s] = projfwd(proj, lat, lon);
    Sen_pos = [x_s(:), y_s(:)];
    N = size(Sen_pos,1);

    paises = shaperead("ne_110m_admin_0_countries.shp", "UseGeoCoords", true);
    colombia = paises(33);

    [xc, yc] = projfwd(proj, colombia.Lat, colombia.Lon);
    polyCol  = build_polyshape(xc, yc);

    %% ================= MALLA (solo dentro de Colombia) =================
    [xmin, xmax] = bounds(polyCol.Vertices(:,1), 'omitnan');
    [ymin, ymax] = bounds(polyCol.Vertices(:,2), 'omitnan');
    x_vec = linspace(xmin, xmax, grid_res);
    y_vec = linspace(ymin, ymax, grid_res);
    [X, Y] = meshgrid(x_vec, y_vec);
    dentro = isinterior(polyCol, X(:), Y(:));
    dentro = reshape(dentro, size(X));

    Err_TDOA = nan(size(X));
    Err_TOA  = nan(size(X));

    %% ================= PROCESAMIENTO =================
    fprintf('Calculando %d puntos dentro de Colombia (x%d MC)... ', nnz(dentro), n_mc);
    t_ini = tic;
    P0 = mean(Sen_pos, 1);                       % semilla inicial

    for i = find(dentro(:))'
        P_real = [X(i), Y(i)];
        d = sqrt(sum((Sen_pos - P_real).^2, 2));
        

        e1 = zeros(n_mc,1); e2 = e1;
        for m = 1:n_mc
            TOA_meas   = d./v_luz + randn(N,1)*sigma_t;
           

            P_tdoa = solve_TDOA_GN(Sen_pos, TOA_meas, v_luz, P0);
            P_toa  = solve_TOA_GN(Sen_pos, TOA_meas, v_luz, P0);
           

            e1(m) = norm(P_real - P_tdoa);
            e2(m) = norm(P_real - P_toa);
            
        end
        Err_TDOA(i) = mean(e1.^2);         % Absoluto
        Err_TOA(i)  = mean(e2.^2);
    end
    fprintf('Hecho (%.1f s).\n', toc(t_ini));

    %% ================= VISUALIZACIÓN =================
    figure('Color','w','Name','Mapas de calor - Red LINET Colombia', ...
           'Position',[60 80 1600 560]);
    tiledlayout(1,3,'TileSpacing','compact','Padding','compact');

    nexttile
    plot_error_map(x_vec, y_vec, Err_TDOA, polyCol, x_s, y_s, ...
        {'TDOA','Diferencia de tiempos de llegada'}, 5000);
    nexttile
    plot_error_map(x_vec, y_vec, Err_TOA, polyCol, x_s, y_s, ...
        {'TOA','Tiempo absoluto (t0 desconocido)'}, 5000);
    

    sgtitle('Error de geolocalización de rayos sobre Colombia - 17 sensores', ...
            'FontSize',14,'FontWeight','bold');

    fprintf('\nError RMS medio dentro de Colombia:\n');
    fprintf('  TDOA : %.1f m\n', mean(Err_TDOA(:),'omitnan'));
    fprintf('  TOA  : %.1f m\n', mean(Err_TOA(:), 'omitnan'));

end

%% ==========================================================================
%  GRAFICADOR: croquis debajo + mapa de calor transparente + contorno encima
%% ==========================================================================
function plot_error_map(x_vec, y_vec, Err, polyCol, x_s, y_s, titulo, cmax)
    hold on
    % 1) Croquis de Colombia (capa inferior)
    plot(polyCol, 'FaceColor',[0.93 0.93 0.93], 'EdgeColor','none', 'FaceAlpha',1);

    % 2) Mapa de calor (NaN = transparente)
    h = imagesc(x_vec, y_vec, Err);
    set(h, 'AlphaData', double(~isnan(Err))*0.9);
    set(gca, 'YDir','normal');
    colormap(gca, turbo);
    clim([0 cmax]);                       % si tu MATLAB es antiguo: caxis([0 cmax])
    cb = colorbar; cb.Label.String = 'Error Absoluto [m]';

    % 3) Contorno encima para bordes limpios
    plot(polyCol, 'FaceColor','none', 'EdgeColor',[0.1 0.1 0.1], 'LineWidth',1.2);

    % 4) Sensores
    plot(x_s, y_s, '^', 'MarkerFaceColor','w', 'MarkerEdgeColor','k', 'MarkerSize',6);

    axis equal tight
    grid on; box on
    xlabel('x (m)'); ylabel('y (m)');
    title(titulo);
    hold off
end

%% ==========================================================================
%  Construye polyshape a partir de contorno con NaN (multipolígono)
%% ==========================================================================
function pol = build_polyshape(xc, yc)
    warning('off','MATLAB:polyshape:repairedBySimplify');
    xc = xc(:)'; yc = yc(:)';
    idxFin = [find(isnan(xc)), numel(xc)+1];
    idxIni = [1, idxFin(1:end-1)+1];
    pol = polyshape();
    for k = 1:numel(idxIni)
        xs = xc(idxIni(k):idxFin(k)-1);
        ys = yc(idxIni(k):idxFin(k)-1);
        if numel(xs) >= 3
            pol = union(pol, polyshape(xs, ys));
        end
    end
end

%% ==========================================================================
%  TDOA - Gauss-Newton (incógnitas x,y). Referencia: sensor 1
%% ==========================================================================
function P = solve_TDOA_GN(S, TOA, c, P0)
    N  = size(S,1);
    dt = TOA - TOA(1);
    idx = 2:N;
    P = P0(:)';
    for it = 1:20
        dist = sqrt(sum((P - S).^2, 2));
        r  = c*dt(idx) - (dist(idx) - dist(1));          % residuo
        u  = (P - S) ./ dist;                            % vectores unitarios
        Jf = u(idx,:) - u(1,:);                          % d(modelo)/dP
        dP = ((Jf'*Jf) + 1e-9*eye(2)) \ (Jf'*r);
        if norm(dP) > 2e5, dP = dP*(2e5/norm(dP)); end   % límite de paso
        P = P + dP';
        if norm(dP) < 1e-3, break; end
    end
end

%% ==========================================================================
%  TOA - Gauss-Newton con incógnitas [x, y, b], donde b = c*t0
%% ==========================================================================
function P = solve_TOA_GN(S, TOA, c, P0)
    N = size(S,1);
    z = c*TOA;
    d0 = sqrt(sum((P0 - S).^2, 2));
    p = [P0(:); mean(z - d0)];                           % [x; y; b]
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

%% ==========================================================================
%  DOA - Intersección por mínimos cuadrados de rectas de acimut
%  Cada sensor: sin(th)*(x-xs) - cos(th)*(y-ys) = 0
%% ==========================================================================
