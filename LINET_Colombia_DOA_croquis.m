function LINET_Colombia_DOA_croquis
% ==========================================================================
% Mapa de calor del error absoluto de geolocalizacion por DOA (acimut)
% sobre Colombia. Red simulada de 17 sensores tipo LINET.
% Proyeccion MAGNA-SIRGAS Origen Nacional (EPSG:9377), con el croquis de
% Colombia debajo del mapa de calor.
% ==========================================================================

    %% ================= PARAMETROS GENERALES =================
    rng(1);
    sigma_theta = deg2rad(3);   % incertidumbre angular [rad]
    grid_res    = 120;            % puntos por eje
    n_mc        = 5;              % repeticiones Monte Carlo por punto
    cmax        = 20000;          % tope de la barra de color [m]

    %% ================= SENSORES =================
    lon = [-77.2703; -75.2686; -71.4211; -76.2824; -74.0727; -75.2166; -73.3448; -74.6447; -73.8388; -75.5545; -70.7452; -72.4610; -73.3188; -72.6169; -74.8266; -75.8665; -76.6724];
    lat = [1.2832; 2.9219; 3.8063; 3.5202; 4.6387; 4.4306; 5.5751; 5.8612; 7.0578; 6.2774; 7.0317; 7.8381; 9.5549; 11.1676; 10.9855; 8.7745; 7.6820];

    %% ================= PROYECCION Y CONTORNO =================
    proj = projcrs(9377);
    [x_s, y_s] = projfwd(proj, lat, lon);
    Sen_pos = [x_s(:), y_s(:)];
    N = size(Sen_pos,1);

    paises = shaperead("ne_110m_admin_0_countries.shp", "UseGeoCoords", true);
    colombia = paises(33);

    [xc, yc] = projfwd(proj, colombia.Lat(:), colombia.Lon(:));
    polyCol  = build_polyshape(xc, yc);

    %% ================= MALLA (solo dentro de Colombia) =================
    [xmin, xmax] = bounds(polyCol.Vertices(:,1), 'omitnan');
    [ymin, ymax] = bounds(polyCol.Vertices(:,2), 'omitnan');
    x_vec = linspace(xmin, xmax, grid_res);
    y_vec = linspace(ymin, ymax, grid_res);
    [X, Y] = meshgrid(x_vec, y_vec);
    dentro = isinterior(polyCol, X(:), Y(:));
    dentro = reshape(dentro, size(X));

    Err_DOA = nan(size(X));

    %% ================= PROCESAMIENTO =================
    fprintf('Calculando %d puntos dentro de Colombia (x%d MC)... ', nnz(dentro), n_mc);
    t_ini = tic;

    for i = find(dentro(:))'
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
    figure('Color','w','Name','Mapa de calor DOA - Red LINET Colombia', ...
           'Position',[100 80 800 700]);
    plot_error_map(x_vec, y_vec, Err_DOA, polyCol, x_s, y_s, ...
        {'DOA','Angulo de llegada (acimut)'}, cmax);

    fprintf('\nError absoluto medio dentro de Colombia:\n');
    fprintf('  DOA : %.1f m\n', mean(Err_DOA(:),'omitnan'));
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
    clim([0 cmax]);                       % MATLAB antiguo: caxis([0 cmax])
    cb = colorbar; cb.Label.String = 'Error absoluto [m]';

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
%  Construye polyshape a partir de contorno con NaN (multipoligono)
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
%  DOA - Gauss-Newton (incognitas x,y) sobre residuos angulares
%  Modelo:    theta_i(P) = atan2(Py - Sy_i, Px - Sx_i)
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

        ang = atan2(dy, dx);
        r = atan2(sin(theta_meas - ang), cos(theta_meas - ang));  % envuelto a [-pi,pi]
        J = [-dy./d2, dx./d2];

        dP = ((J'*J) + 1e-18*eye(2)) \ (J'*r);
        if norm(dP) > 2e5, dP = dP*(2e5/norm(dP)); end   % limite de paso
        P = P + dP';
        if norm(dP) < 1e-3, break; end
    end
end