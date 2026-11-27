function S = mapaErrorColombia(opts)
% MAPAERRORCOLOMBIA  Simula rayos en posiciones CONOCIDAS dentro de Colombia,
% los localiza con el detector detectarRayosGLM sobre una imagen sintetica de la
% region (1 pixel = opts.pixelM metros, por defecto 4 km x 4 km) y grafica el
% mapa de calor del error absoluto [m] con el estilo de plot_error_map
% (croquis abajo + mapa transparente + contorno).
%
% El flash detectado (centroide ponderado, en pixeles) se convierte de nuevo a
% metros y se compara con la posicion verdadera.
%
% USO
%   S = mapaErrorColombia();
%   S = mapaErrorColombia(struct('pixelM',4000,'nRep',2));
%   S = mapaErrorColombia(struct('contorno',[xc(:) yc(:)]));
%
% Genera 2 figuras: (1) mapa de error y (2) rayos simulados vs. estimados.
%
% IMPORTANTE
%   * El contorno por defecto sale del shapefile Natural Earth (opts.shapefile),
%     proyectado con projcrs(9377) (MAGNA-SIRGAS Origen Nacional, metros).
%     Para usar otro, pasa opts.contorno = [xc yc] en metros (con NaN entre
%     poligonos).
%   * Coordenadas en metros segun EPSG:9377.
%
% OPCIONES (struct, todas opcionales)
%   contorno [] (Nx2 m)   shapefile 'ne_110m_admin_0_countries.shp'
%   pixelM 4000      tamano de pixel de la imagen simulada [m]
%   sigmaKm 12       radio de la mancha de luz del rayo [km]
%   espaciadoKm 16   separacion de la grilla regular de rayos [km]; cada celda
%                    de la grilla es una celda del mapa de calor
%   nRep 1           repeticiones por posicion (ruido de imagen independiente)
%   rMaxKm 100       distancia maxima para emparejar flash detectado con verdad
%   A 0.35, fondo 0.10, ruido 0.03   amplitud, fondo y ruido de la imagen (0-1)
%   N 40, kIni 15, nFlash 10         frames totales, inicio y duracion del flash
%   det              struct con opciones extra para detectarRayosGLM
%   cmax []          limite superior de color [m] (vacio = percentil 95)
%   amplificar 20    factor que exagera los segmentos de error en la figura 2
%   semilla 1, csv ''  (csv = prefijo; guarda prefijo_imagen.csv)
%
% Requiere: Image Processing Toolbox, Mapping Toolbox y detectarRayosGLM.m (version con entrada
% numerica y opciones graficar/verbose). polyshape requiere R2017b o superior.

if nargin < 1, opts = struct(); end
s = opcionesCol(opts);
rng(s.semilla);

%% 1. Contorno y dominio ------------------------------------------------------
if isempty(s.contorno)
    [xc, yc] = contornoColombiaShp(s.shapefile);
else
    xc = s.contorno(:,1);  yc = s.contorno(:,2);
end
polyCol = build_polyshape(xc, yc);

[xl, yl] = boundingbox(polyCol);
mrg = 5 * s.pixelM;                                  % margen alrededor del pais
dom.xmin = xl(1) - mrg;  dom.xmax = xl(2) + mrg;
dom.ymin = yl(1) - mrg;  dom.ymax = yl(2) + mrg;
dom.W = ceil((dom.xmax - dom.xmin) / s.pixelM) + 1;
dom.H = ceil((dom.ymax - dom.ymin) / s.pixelM) + 1;

%% 2. Posiciones verdaderas de los rayos: GRILLA REGULAR dentro del contorno ---
% Un rayo en el centro de cada celda de espaciadoKm x espaciadoKm que caiga
% dentro del pais. Cada celda de la grilla es tambien una celda del mapa de calor.
esp = s.espaciadoKm * 1000;
x_vec = (xl(1) + esp/2) : esp : xl(2);
y_vec = (yl(1) + esp/2) : esp : yl(2);
nx = numel(x_vec);  ny = numel(y_vec);
[Xg, Yg] = meshgrid(x_vec, y_vec);
in = isinterior(polyCol, Xg(:), Yg(:));
[IY0, IX0] = ind2sub([ny nx], find(in));
P0 = [x_vec(IX0(:)).' , y_vec(IY0(:)).'];
% nRep > 1 repite cada posicion con ruido de imagen independiente
P  = repmat(P0,  s.nRep, 1);
IX = repmat(IX0(:), s.nRep, 1);
IY = repmat(IY0(:), s.nRep, 1);
fprintf('Grilla de %g km: %d x %d celdas, %d rayos dentro del contorno (x%d rep.)\n', ...
    s.espaciadoKm, nx, ny, size(P0,1), s.nRep);

%% 3. Estimacion de la posicion por imagen -----------------------------------
[xd, yd] = estimarPorImagen(P, dom, s);
nombre = sprintf('Detector por imagen (tipo GLM, pixel %g km)', s.pixelM/1000);

%% 4. Errores absolutos, resumen y mapa --------------------------------------
err = hypot(xd - P(:,1), yd - P(:,2));                   % error absoluto [m]
Err = accumarray([IY IX], err, [ny nx], @(a) mean(a, 'omitnan'), NaN);
resumenErr(nombre, err);

v = Err(~isnan(Err));
if isempty(s.cmax), cmax = percentil(v, 0.95); else, cmax = s.cmax; end
if isempty(cmax) || ~(cmax > 0), cmax = 1; end

figure('Name', 'Mapa de calor del error absoluto', 'Color', 'w', ...
       'Position', [80 80 560 640]);
plot_error_map(x_vec, y_vec, Err, polyCol, nombre, cmax);

T = table(P(:,1), P(:,2), xd, yd, err, 'VariableNames', ...
          {'x_verdad_m','y_verdad_m','x_est_m','y_est_m','err_abs_m'});
S = struct('polyCol', polyCol, 'x_vec', x_vec, 'y_vec', y_vec, ...
           'imagen', struct('tabla', T, 'Err', Err));
if ~isempty(s.csv), writetable(T, [s.csv '_imagen.csv']); end

%% 5. Figura 2: rayos simulados vs. estimados --------------------------------
figure('Name', 'Rayos simulados', 'Color', 'w', 'Position', [100 100 560 640]);
plot_rayos(polyCol, P, xd, yd, nombre, s.amplificar);
end

%% ========================= OPCIONES Y DATOS ==============================

function s = opcionesCol(opts)
s = struct('contorno',[], 'shapefile','ne_110m_admin_0_countries.shp', 'pixelM',4000, ...
    'sigmaKm',12, 'espaciadoKm',16, 'nRep',1, 'rMaxKm',100, ...
    'A',0.35, 'fondo',0.10, 'ruido',0.03, 'N',40, 'kIni',15, 'nFlash',10, ...
    'det',struct(), 'cmax',[], 'amplificar',20, 'semilla',1, 'csv','');
f = fieldnames(opts);
for i = 1:numel(f), s.(f{i}) = opts.(f{i}); end
end

function [xc, yc] = contornoColombiaShp(archivo)
% Contorno de Colombia desde Natural Earth, proyectado a MAGNA-SIRGAS Origen
% Nacional (EPSG:9377, metros). Requiere Mapping Toolbox.
proj = projcrs(9377);
paises = shaperead(archivo, "UseGeoCoords", true);
idx = [];
for campo = ["NAME","ADMIN","NAME_LONG","SOVEREIGNT"]
    if isfield(paises, campo)
        idx = find(strcmpi({paises.(campo)}, 'Colombia'), 1);
        if ~isempty(idx), break; end
    end
end
if isempty(idx), idx = 33; end          % respaldo: indice usado antes
colombia = paises(idx);
[xc, yc] = projfwd(proj, colombia.Lat(:), colombia.Lon(:));
end

%% ========================= ESTIMADOR =====================================

function [xd, yd] = estimarPorImagen(P, dom, s)
% Cada rayo se dibuja como una mancha gaussiana en una imagen de la region, con
% ruido, y se procesa con detectarRayosGLM. Pixel (1,1) = esquina NO del dominio.
n = size(P,1);  xd = nan(n,1);  yd = nan(n,1);
[X, Y] = meshgrid(1:dom.W, 1:dom.H);
sig = s.sigmaKm * 1000 / s.pixelM;                    % radio en pixeles
rMaxPx = s.rMaxKm * 1000 / s.pixelM;
% maxDistPx se escala con el tamano de pixel (equivale a ~320 km)
dopt = struct('graficar', false, 'verbose', false, 'minPixGrupo', 3, ...
              'minGruposFlash', 2, 'maxDistPx', round(320000 / s.pixelM));
f = fieldnames(s.det);
for i = 1:numel(f), dopt.(f{i}) = s.det.(f{i}); end
fprintf('Simulando %d rayos por imagen (%dx%d px)...\n', n, dom.W, dom.H);
for k = 1:n
    px = (P(k,1) - dom.xmin) / s.pixelM + 1;          % metros -> pixeles
    py = (dom.ymax - P(k,2)) / s.pixelM + 1;          % norte arriba
    F = s.fondo + s.ruido * randn(dom.H, dom.W, s.N, 'single');
    B = s.A * exp(-((X - px).^2 + (Y - py).^2) / (2 * sig^2));
    for q = 1:s.nFlash
        kk = s.kIni + q - 1;
        F(:,:,kk) = F(:,:,kk) + sin(pi * q / (s.nFlash + 1)) * B;
    end
    R = detectarRayosGLM(F, dopt);
    if height(R.flashes) > 0
        d = hypot(R.flashes.x_px - px, R.flashes.y_px - py);
        [dm, j] = min(d);
        if dm <= rMaxPx
            xd(k) = dom.xmin + (R.flashes.x_px(j) - 1) * s.pixelM;   % pixeles -> metros
            yd(k) = dom.ymax - (R.flashes.y_px(j) - 1) * s.pixelM;
        end
    end
    if mod(k, 50) == 0, fprintf('  %d/%d\n', k, n); end
end
end

%% ========================= ESTADISTICAS ==================================

function resumenErr(nombre, e)
ok = ~isnan(e);
fprintf('\n=== %s ===\n', nombre);
fprintf('Rayos: %d | localizados: %d (%.1f %%)\n', numel(e), nnz(ok), 100 * mean(ok));
if ~any(ok), disp('Ninguno localizado.'); return; end
v = sort(e(ok));
fprintf('MAE = %.1f m | RMSE = %.1f m | mediana = %.1f m | P95 = %.1f m | max = %.1f m\n', ...
    mean(v), sqrt(mean(v.^2)), median(v), percentil(v, 0.95), v(end));
end

function p = percentil(v, q)
v = sort(v(:));
if isempty(v), p = []; return; end
p = v(max(1, ceil(q * numel(v))));
end

%% ==========================================================================
%  GRAFICADOR: croquis debajo + mapa de calor transparente + contorno encima
%% ==========================================================================
function plot_error_map(x_vec, y_vec, Err, polyCol, titulo, cmax)
    hold on
    % 1) Croquis de Colombia (capa inferior)
    plot(polyCol, 'FaceColor',[0.93 0.93 0.93], 'EdgeColor','none', 'FaceAlpha',1);
    % 2) Mapa de calor (NaN = transparente)
    h = imagesc(x_vec, y_vec, Err);
    set(h, 'AlphaData', double(~isnan(Err))*0.9);
    set(gca, 'YDir','normal');
    colormap(gca, turbo);
    try, clim([0 cmax]); catch, caxis([0 cmax]); end   %#ok<CAXIS>
    cb = colorbar; cb.Label.String = 'Error absoluto [m]';
    % 3) Contorno encima para bordes limpios
    plot(polyCol, 'FaceColor','none', 'EdgeColor',[0.1 0.1 0.1], 'LineWidth',1.2);
    axis equal tight
    grid on; box on
    xlabel('x (m)'); ylabel('y (m)');
    title(titulo);
    hold off
end

%% Rayos simulados (verdad) y posicion estimada, con el mismo croquis
function plot_rayos(polyCol, P, xd, yd, titulo, amp)
    hold on
    plot(polyCol, 'FaceColor',[0.93 0.93 0.93], 'EdgeColor',[0.1 0.1 0.1], ...
         'LineWidth',1.2, 'FaceAlpha',1);
    xa = P(:,1) + amp * (xd - P(:,1));            % error exagerado para verlo
    ya = P(:,2) + amp * (yd - P(:,2));
    L = [P(:,1), xa, nan(size(xa))]';
    M = [P(:,2), ya, nan(size(ya))]';
    plot(L(:), M(:), '-', 'Color',[0.85 0.2 0.2]);
    h1 = plot(P(:,1), P(:,2), '.', 'Color',[0.1 0.3 0.8], 'MarkerSize',9);
    h2 = plot(xa, ya, 'x', 'Color',[0.85 0.2 0.2], 'MarkerSize',4);
    legend([h1 h2], {'Rayo simulado (verdad)', 'Posicion estimada'}, 'Location','best');
    axis equal tight; grid on; box on
    xlabel('x (m)'); ylabel('y (m)');
    title(sprintf('%s (error ampliado x%g)', titulo, amp));
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