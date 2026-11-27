function R = detectarRayosGLM(entrada, opts)
% DETECTARRAYOSGLM  Detecta rayos con la logica de la guia GLM (GOES-R):
%   Evento  -> pixel que supera el umbral respecto a un fondo que se actualiza
%   Grupo   -> eventos simultaneos en pixeles vecinos (conectividad 8)
%   Flash   -> grupos consecutivos separados < maxGapMs y < maxDistPx
%   Productos: FED (densidad de extension de flashes), AFA (area media) y
%   TOE (energia optica total), igual que los "gridded products" del PDF.
%
% USO
%   R = detectarRayosGLM('demo');                    % prueba con datos sinteticos
%   R = detectarRayosGLM('foto.jpg');                % UNA imagen (modo espacial)
%   R = detectarRayosGLM('video.mp4');               % video (modo temporal)
%   R = detectarRayosGLM('carpeta_frames');          % secuencia de imagenes
%   R = detectarRayosGLM('video.mp4', struct('K',4,'maxDistPx',60));
%   R = detectarRayosGLM(F, opts);                   % F = arreglo numerico HxWxN (pruebas)
%
% OPCIONES (struct, todas opcionales)
%   K                 umbral = K * sigma_ruido            (def. 5)
%   umbralMin         umbral minimo absoluto, imagen 0-1  (def. 0.05)
%   alpha             tasa de actualizacion del fondo     (def. 0.05)
%   nInicioBg         frames iniciales para el fondo      (def. 10)
%   maxFraccionEventos si > fraccion de pixeles se activa, se descarta el
%                     frame (perturbacion global/plataforma) (def. 0.2)
%   maxPersistencia   pixel activo en > esta fraccion de frames = falso
%                     (pixel caliente, reflejo, glint)    (def. 0.5)
%   minPixGrupo       pixeles minimos por grupo           (def. 1)
%   minGruposFlash    grupos minimos por flash            (def. 1)
%   dtMs              tiempo entre frames en ms (GLM = 2) (def. 2)
%   maxGapMs          separacion temporal maxima (GLM=330)(def. 330)
%   maxDistPx         separacion espacial maxima en pixeles (def. 40)
%                     (GLM: 16.5 km / 8 km por pixel ~ 2 px)
%   pixKm             km por pixel (GLM ~8; def. 4 = areas en px^2)
%   escala            reduce la imagen (0.5 = mitad)      (def. 1)
%   maxFrames         tope de frames a cargar             (def. 3000)
%   radioFondo        radio top-hat, modo foto unica      (def. 15)
%   solidezMax, excentricidadMin, minPixLargo  filtro de forma, foto unica
%   csv               ruta de un .csv para guardar los flashes (def. '')
%   graficar          true/false: mostrar figuras          (def. true)
%   verbose           true/false: imprimir mensajes        (def. true)
%
% Requiere: Image Processing Toolbox.

if nargin < 1, entrada = 'demo'; end
if nargin < 2, opts = struct(); end
o = opcionesPorDefecto(opts);

%% 1. Cargar datos --------------------------------------------------------
if isnumeric(entrada)                      % arreglo HxWxN ya cargado en memoria
    F = single(entrada); tms = (0:size(F,3)-1) * o.dtMs;
elseif strcmpi(entrada, 'demo')
    [F, tms] = datosDemo(o);
else
    [F, tms, o] = cargarFrames(entrada, o);
end
[H, W, N] = size(F);
esSecuencia = N > 1;
imprimir(o, 'Cargado: %d frame(s) de %dx%d px. Modo: %s\n', N, W, H, ...
    ternario(esSecuencia, 'temporal (fondo que se actualiza)', 'foto unica (fondo espacial)'));

%% 2. EVENTOS: pixeles sobre el umbral respecto al fondo -------------------
if esSecuencia
    [EV, D] = eventosTemporal(F, o);
else
    [EV, D] = eventosFotoUnica(F, o);
end
imprimir(o, 'Eventos (pixel-frame): %d\n', nnz(EV));

%% 3. GRUPOS: eventos vecinos en el mismo frame ---------------------------
gT = []; gX = []; gY = []; gE = []; gN = []; gIdx = {}; gVal = {};
for k = 1:N
    if ~any(any(EV(:,:,k))), continue; end
    cc = bwconncomp(EV(:,:,k), 8);
    Dk = D(:,:,k);
    for j = 1:cc.NumObjects
        p = cc.PixelIdxList{j};
        if numel(p) < o.minPixGrupo, continue; end
        w = double(Dk(p));
        [r, c] = ind2sub([H W], p);
        gT(end+1,1) = tms(k);                       %#ok<AGROW>
        gX(end+1,1) = sum(c .* w) / sum(w);         %#ok<AGROW> centroide ponderado
        gY(end+1,1) = sum(r .* w) / sum(w);         %#ok<AGROW> por radiancia
        gE(end+1,1) = sum(w);                       %#ok<AGROW> energia del grupo
        gN(end+1,1) = numel(p);                     %#ok<AGROW>
        gIdx{end+1,1} = p;                          %#ok<AGROW>
        gVal{end+1,1} = w;                          %#ok<AGROW>
    end
end
nG = numel(gT);
imprimir(o, 'Grupos: %d\n', nG);

R = struct('flashes', table(), 'FED', zeros(H,W), 'AFA', nan(H,W), ...
           'TOE', zeros(H,W), 'EV', EV);
if nG == 0
    imprimir(o, 'No se detectaron grupos. Prueba bajar K o umbralMin.\n');
    return
end

%% 4. FLASHES: grupos consecutivos cerca en tiempo y espacio ---------------
flashId = zeros(nG,1); nF = 0; lastT = zeros(0,1);
for g = 1:nG
    mejor = 0; dMin = inf;
    for f = 1:nF
        if gT(g) - lastT(f) > o.maxGapMs, continue; end
        m = find(flashId(1:g-1) == f & (gT(g) - gT(1:g-1)) <= o.maxGapMs);
        d = min(hypot(gX(m) - gX(g), gY(m) - gY(g)));
        if d < o.maxDistPx && d < dMin, mejor = f; dMin = d; end
    end
    if mejor == 0
        nF = nF + 1; mejor = nF; lastT(nF,1) = gT(g);
    else
        lastT(mejor) = gT(g);
    end
    flashId(g) = mejor;
end

ids = find(arrayfun(@(f) nnz(flashId == f) >= o.minGruposFlash, 1:nF));
nK = numel(ids);
imprimir(o, 'Flashes: %d\n', nK);
if nK == 0, return; end

%% 5. Estadisticas por flash y productos en rejilla (FED, AFA, TOE) --------
FED = zeros(H,W); AFAs = zeros(H,W); TOEv = zeros(H*W,1);
fx = zeros(nK,1); fy = fx; t0 = fx; dur = fx; ng = fx; ne = fx; ar = fx; en = fx;
for i = 1:nK
    m = find(flashId == ids(i));
    w = gE(m);
    fx(i) = sum(gX(m) .* w) / sum(w);               % ubicacion del flash =
    fy(i) = sum(gY(m) .* w) / sum(w);               % centroide ponderado
    idxAll = vertcat(gIdx{m});  valAll = vertcat(gVal{m});
    [u, ~, ic] = unique(idxAll);                    % huella del flash
    ar(i) = numel(u) * o.pixKm^2;
    t0(i) = min(gT(m));
    dur(i) = max(gT(m)) - t0(i) + o.dtMs;
    ng(i) = numel(m);  ne(i) = sum(gN(m));  en(i) = sum(w);
    FED(u)  = FED(u) + 1;
    AFAs(u) = AFAs(u) + ar(i);
    TOEv(u) = TOEv(u) + accumarray(ic, double(valAll));
end
AFA = AFAs ./ FED;  AFA(FED == 0) = NaN;
TOE = reshape(TOEv, H, W);

R.flashes = table((1:nK)', t0, dur, fx, fy, ng, ne, ar, en, 'VariableNames', ...
    {'ID','t0_ms','duracion_ms','x_px','y_px','n_grupos','n_eventos','area','energia'});
R.FED = FED; R.AFA = AFA; R.TOE = TOE;
if o.verbose, disp(R.flashes); end
if ~isempty(o.csv), writetable(R.flashes, o.csv); imprimir(o, 'Guardado: %s\n', o.csv); end

%% 6. Visualizacion --------------------------------------------------------
if o.graficar
ref = mat2gray(max(F, [], 3));
refRGB = repmat(ref, [1 1 3]);
gKeep = ismember(flashId, ids);

figure('Name', 'Deteccion de rayos (estilo GLM)', 'Color', 'w');
subplot(2,2,1); imshow(refRGB); hold on;
plot(gX(gKeep), gY(gKeep), 'w.', 'MarkerSize', 6);
plot(fx, fy, 'gx', 'MarkerSize', 12, 'LineWidth', 2);
title(sprintf('Grupos (blanco) y flashes (X verde): %d', nK));
subplot(2,2,2); superponer(refRGB, FED, 'FED (flashes por pixel)', false);
subplot(2,2,3); superponer(refRGB, AFA, 'AFA (area media del flash)', false);
subplot(2,2,4); superponer(refRGB, TOE, 'TOE (energia optica, log)', true);

if esSecuencia
    figure('Name', 'Serie temporal', 'Color', 'w');
    ev_t = squeeze(sum(sum(EV, 1), 2));
    bar(tms, ev_t, 1); hold on;
    for i = 1:nK
        xline(t0(i), 'g'); 
    end
    xlabel('tiempo (ms)'); ylabel('eventos por frame');
    title('Emisiones opticas discretas (lineas verdes = inicio de flash)');
end
end   % cierra if o.graficar
end

%% ========================= FUNCIONES LOCALES ============================

function o = opcionesPorDefecto(opts)
o = struct('K',5, 'umbralMin',0.05, 'alpha',0.05, 'nInicioBg',10, ...
    'maxFraccionEventos',0.2, 'maxPersistencia',0.5, 'minPixGrupo',1, ...
    'minGruposFlash',1, 'dtMs',2, 'maxGapMs',330, 'maxDistPx',40, ...
    'pixKm',1, 'escala',1, 'maxFrames',3000, 'radioFondo',15, ...
    'solidezMax',0.6, 'excentricidadMin',0.97, 'minPixLargo',15, 'csv','', ...
    'graficar',true, 'verbose',true);
f = fieldnames(opts);
for i = 1:numel(f), o.(f{i}) = opts.(f{i}); end
end

function [EV, D] = eventosTemporal(F, o)
% Fondo que se actualiza continuamente; el fondo NO se actualiza donde hay
% evento, para que el rayo no "contamine" el fondo.
[H, W, N] = size(F);
bg = median(F(:,:,1:min(o.nInicioBg, N)), 3);
EV = false(H, W, N); D = zeros(H, W, N, 'single');
for k = 1:N
    I = F(:,:,k);
    d = I - bg;
    sigma = max(1.4826 * median(abs(d(:) - median(d(:)))), 1e-3);  % ruido robusto (MAD)
    ev = d > max(o.K * sigma, o.umbralMin);
    if nnz(ev) / numel(ev) > o.maxFraccionEventos   % cambio global: no es rayo
        ev(:) = false; bg = I; 
    else
        bg(~ev) = bg(~ev) + o.alpha * (I(~ev) - bg(~ev));
    end
    EV(:,:,k) = ev; D(:,:,k) = d .* ev;
end
% Filtro de falsos eventos persistentes (pixel caliente, glint, reflejo fijo)
if N >= 20
    frec = mean(EV, 3);
    EV(repmat(frec > o.maxPersistencia, [1 1 N])) = false;
    D = D .* EV;
end
end

function [EV, D] = eventosFotoUnica(I, o)
% Sin tiempo no hay "cambio": el fondo se estima espacialmente (top-hat) y
% se filtran formas (el rayo es delgado/ramificado; estrellas y luces son
% manchas compactas).
th = imtophat(I, strel('disk', o.radioFondo));
sigma = max(1.4826 * median(abs(th(:) - median(th(:)))), 1e-3);
ev = th > max(o.K * sigma, o.umbralMin);
cc = bwconncomp(ev, 8);
EV = false(size(I)); 
if cc.NumObjects > 0
    st = regionprops(cc, 'Solidity', 'Eccentricity', 'Area');
    keep = ([st.Solidity] < o.solidezMax) | ...
           ([st.Eccentricity] > o.excentricidadMin & [st.Area] >= o.minPixLargo);
    EV(vertcat(cc.PixelIdxList{keep})) = true;
end
D = th .* EV;
end

function [F, t, o] = cargarFrames(entrada, o)
if isfolder(entrada)
    d = dir(entrada); d = d(~[d.isdir]);
    [~, ~, e] = cellfun(@fileparts, {d.name}, 'UniformOutput', false);
    d = d(ismember(lower(e), {'.jpg','.jpeg','.png','.tif','.tiff','.bmp'}));
    [~, ord] = sort({d.name}); d = d(ord);
    N = min(numel(d), o.maxFrames); c = cell(1, N);
    for k = 1:N, c{k} = leerGris(fullfile(d(k).folder, d(k).name), o.escala); end
    F = cat(3, c{:}); t = (0:N-1) * o.dtMs;
else
    [~, ~, ext] = fileparts(entrada);
    if any(strcmpi(ext, {'.mp4','.avi','.mov','.mkv','.m4v'}))
        v = VideoReader(entrada);
        o.dtMs = 1000 / v.FrameRate;          % ms reales entre frames
        c = {};
        while hasFrame(v) && numel(c) < o.maxFrames
            I = readFrame(v);
            c{end+1} = leerGrisArray(I, o.escala); %#ok<AGROW>
        end
        F = cat(3, c{:}); t = (0:numel(c)-1) * o.dtMs;
    else
        F = leerGris(entrada, o.escala); t = 0;
    end
end
end

function I = leerGris(archivo, esc)
I = leerGrisArray(imread(archivo), esc);
end

function I = leerGrisArray(I, esc)
if size(I,3) == 3, I = rgb2gray(I); end
I = im2single(I);
if esc ~= 1, I = imresize(I, esc); end
end

function superponer(refRGB, M, titulo, usarLog)
imshow(refRGB); hold on;
if usarLog, M = log10(M + 1); end
h = imagesc(M);
set(h, 'AlphaData', 0.75 * double(~isnan(M) & M > 0));
colormap(gca, parula); colorbar; title(titulo);
end

function imprimir(o, varargin)
if o.verbose, fprintf(varargin{:}); end
end

function s = ternario(c, a, b)
if c, s = a; else, s = b; end
end

%% ---------------------- DATOS SINTETICOS DE PRUEBA ----------------------
function [F, t] = datosDemo(o)
rng(1); H = 240; W = 320; N = 300;
[X, Y] = meshgrid(1:W, 1:H);
nube = 0.15 * exp(-((X-160).^2 + (Y-100).^2) / (2*60^2));
F = 0.10 + nube + 0.02 * randn(H, W, N, 'single');
% Flash A: canal que crece cerca de (80,60), con una pausa de "oscuridad"
for k = 50:70
    if k >= 58 && k <= 60, continue; end
    n = (k - 49) * 3;
    F = dibujarCanal(F, k, 80, 60, 80 + n, 60 + round(n * 0.6), 0.4);
end
% Flash B: lejos de A (mas de maxDistPx) -> debe ser otro flash
for k = 150:175
    n = (k - 149) * 3;
    F = dibujarCanal(F, k, 240, 170, 240 - n, 170 - round(n * 0.4), 0.45);
end
% Falso evento: pixel caliente que se enciende y queda fijo (filtro de persistencia)
F(30:31, 30:31, 100:N) = 0.7;
t = (0:N-1) * o.dtMs;
end

function F = dibujarCanal(F, k, x0, y0, x1, y1, br)
n = max(abs(x1-x0), abs(y1-y0)) + 1;
x = min(max(round(linspace(x0, x1, n)), 2), size(F,2)-1);
y = min(max(round(linspace(y0, y1, n)), 2), size(F,1)-1);
Ik = F(:,:,k);
for dx = -1:1
    for dy = 0:1
        ind = sub2ind(size(Ik), y + dy, x + dx);
        Ik(ind) = Ik(ind) + br;
    end
end
F(:,:,k) = Ik;
end
