clear; clc; close all;

rng(1);

rootDir = fileparts(mfilename('fullpath'));

addpath(fullfile(rootDir, 'polinbe'));
addpath(fullfile(rootDir, 'polinbe', 'analytical_integrals'));

nCasos = 100;
nMallasMostrar = 4;

Ta = 50.0;
hConv = 60.0;
qBase = -23292.0;

W = 0.028;
Hbase = 0.005;
HtotalMax = 0.025;

nAletas = 10;
espesor = 0.001;
separacion = 0.002;

Le = 0.001;

hAletaMin = 0.005;
hAletaMax = HtotalMax - Hbase;

materiales = {
    'Aluminio', 200.0;
    'Acero',     50.0
};

carpetaCasos = fullfile(rootDir, 'casos_10_aletas');

if ~exist(carpetaCasos, 'dir')
    mkdir(carpetaCasos);
end

carpetaFiguras = fullfile(rootDir, 'figuras');

if ~exist(carpetaFiguras, 'dir')
    mkdir(carpetaFiguras);
end

todosResultados = struct();

for imat = 1:size(materiales, 1)

    nombreMaterial = materiales{imat, 1};
    k = materiales{imat, 2};

    Lconv = zeros(nCasos, 1);
    TbaseProm = zeros(nCasos, 1);
    alturasCasos = zeros(nCasos, nAletas);

    fprintf('\nMaterial: %s, k = %.1f W/(m K)\n', nombreMaterial, k);

    for icase = 1:nCasos

        hAletas = hAletaMin + (hAletaMax - hAletaMin) * rand(1, nAletas);
        alturasCasos(icase, :) = hAletas;

        infile = fullfile(carpetaCasos, sprintf('modelo_%s_%04d.txt', nombreMaterial, icase));
        outfile = fullfile(carpetaCasos, sprintf('salida_%s_%04d.txt', nombreMaterial, icase));

        [nodes, elems] = generar_modelo_10_aletas( ...
            infile, Ta, k, hConv, qBase, W, Hbase, ...
            nAletas, espesor, separacion, hAletas, Le);

        Lconv(icase) = calcular_longitud_convectiva(nodes, elems);

        figurasAntes = findall(0, 'Type', 'figure');

        if imat == 1 && icase <= nMallasMostrar
            set(0, 'DefaultFigureVisible', 'on');
        else
            set(0, 'DefaultFigureVisible', 'off');
        end

        polinbe_cc(infile, outfile);

        set(0, 'DefaultFigureVisible', 'on');

        figurasDespues = findall(0, 'Type', 'figure');
        figurasNuevas = setdiff(figurasDespues, figurasAntes);

        if ~(imat == 1 && icase <= nMallasMostrar)
            close(figurasNuevas);
        end

        datos = leer_salida_polinbe(outfile);

        TbaseProm(icase) = calcular_temperatura_media_base(datos);

        fprintf('%s caso %4d/%4d | Lconv = %.5f m | Tbase = %.3f ºC\n', ...
            nombreMaterial, icase, nCasos, Lconv(icase), TbaseProm(icase));
    end

    todosResultados(imat).material = nombreMaterial;
    todosResultados(imat).k = k;
    todosResultados(imat).Lconv = Lconv;
    todosResultados(imat).TbaseProm = TbaseProm;
    todosResultados(imat).alturas = alturasCasos;

    tablaResultados = array2table(alturasCasos);

    for i = 1:nAletas
        tablaResultados.Properties.VariableNames{i} = sprintf('h_aleta_%02d_m', i);
    end

    tablaResultados.Lconv_m = Lconv;
    tablaResultados.TbaseProm_C = TbaseProm;

    writetable(tablaResultados, fullfile(carpetaCasos, sprintf('resultados_%s.csv', nombreMaterial)));

    analizar_y_graficar_material(nombreMaterial, Lconv, TbaseProm, Ta, rootDir, carpetaFiguras);
end

save(fullfile(carpetaCasos, 'resultados_completos.mat'), 'todosResultados');


function [nodes, elems] = generar_modelo_10_aletas( ...
    filename, Ta, k, hConv, qBase, W, Hbase, ...
    nAletas, espesor, separacion, hAletas, Le)

    contornos = [
        1, 0.0;
        2, hConv;
        1, qBase
    ];

    BID_ADIABATICO = 1;
    BID_CONVECCION = 2;
    BID_BASE = 3;

    segmentos = {};

    segmentos{end + 1} = struct( ...
        'p0', [0.0, 0.0], ...
        'p1', [W, 0.0], ...
        'bid', BID_BASE);

    xL = zeros(1, nAletas);
    xR = zeros(1, nAletas);
    yTop = zeros(1, nAletas);

    for i = 1:nAletas
        xL(i) = (i - 1) * (espesor + separacion);
        xR(i) = xL(i) + espesor;
        yTop(i) = Hbase + hAletas(i);
    end

    segmentos{end + 1} = struct( ...
        'p0', [W, 0.0], ...
        'p1', [W, yTop(end)], ...
        'bid', BID_ADIABATICO);

    segmentos{end + 1} = struct( ...
        'p0', [xR(end), yTop(end)], ...
        'p1', [xL(end), yTop(end)], ...
        'bid', BID_CONVECCION);

    segmentos{end + 1} = struct( ...
        'p0', [xL(end), yTop(end)], ...
        'p1', [xL(end), Hbase], ...
        'bid', BID_CONVECCION);

    for i = nAletas - 1:-1:2

        segmentos{end + 1} = struct( ...
            'p0', [xL(i + 1), Hbase], ...
            'p1', [xR(i), Hbase], ...
            'bid', BID_CONVECCION);

        segmentos{end + 1} = struct( ...
            'p0', [xR(i), Hbase], ...
            'p1', [xR(i), yTop(i)], ...
            'bid', BID_CONVECCION);

        segmentos{end + 1} = struct( ...
            'p0', [xR(i), yTop(i)], ...
            'p1', [xL(i), yTop(i)], ...
            'bid', BID_CONVECCION);

        segmentos{end + 1} = struct( ...
            'p0', [xL(i), yTop(i)], ...
            'p1', [xL(i), Hbase], ...
            'bid', BID_CONVECCION);
    end

    segmentos{end + 1} = struct( ...
        'p0', [xL(2), Hbase], ...
        'p1', [xR(1), Hbase], ...
        'bid', BID_CONVECCION);

    segmentos{end + 1} = struct( ...
        'p0', [xR(1), Hbase], ...
        'p1', [xR(1), yTop(1)], ...
        'bid', BID_CONVECCION);

    segmentos{end + 1} = struct( ...
        'p0', [xR(1), yTop(1)], ...
        'p1', [xL(1), yTop(1)], ...
        'bid', BID_CONVECCION);

    segmentos{end + 1} = struct( ...
        'p0', [xL(1), yTop(1)], ...
        'p1', [0.0, 0.0], ...
        'bid', BID_ADIABATICO);

    nodes = [];
    elems = [];

    for iseg = 1:numel(segmentos)

        p0 = segmentos{iseg}.p0;
        p1 = segmentos{iseg}.p1;
        bid = segmentos{iseg}.bid;

        pts = dividir_segmento(p0, p1, Le);

        start = size(nodes, 1) + 1;

        nodes = [nodes; pts];

        for j = 1:size(pts, 1) - 1
            elems = [elems; start + j - 1, start + j, bid];
        end
    end

    internal = generar_puntos_internos(W, Hbase);

    fid = fopen(filename, 'w');

    if fid == -1
        error('No se pudo crear el fichero %s', filename);
    end

    fprintf(fid, '%.12g\n', Ta);
    fprintf(fid, '%.12g\n', k);

    fprintf(fid, '%d\n', size(contornos, 1));
    for i = 1:size(contornos, 1)
        fprintf(fid, '%d %.12g\n', contornos(i, 1), contornos(i, 2));
    end

    fprintf(fid, '%d\n', size(nodes, 1));
    for i = 1:size(nodes, 1)
        fprintf(fid, '%.12g %.12g\n', nodes(i, 1), nodes(i, 2));
    end

    fprintf(fid, '%d\n', size(elems, 1));
    for i = 1:size(elems, 1)
        fprintf(fid, '%d %d %d\n', elems(i, 1), elems(i, 2), elems(i, 3));
    end

    fprintf(fid, '%d\n', size(internal, 1));
    for i = 1:size(internal, 1)
        fprintf(fid, '%.12g %.12g\n', internal(i, 1), internal(i, 2));
    end

    fclose(fid);
end


function pts = dividir_segmento(p0, p1, Le)

    dx = p1(1) - p0(1);
    dy = p1(2) - p0(2);

    longitud = hypot(dx, dy);
    n = max(1, round(longitud / Le));

    pts = zeros(n + 1, 2);

    for i = 0:n
        pts(i + 1, :) = p0 + (p1 - p0) * i / n;
    end
end


function internal = generar_puntos_internos(W, Hbase)

    xMid = W / 2;
    yMin = 0.0005;
    yMax = Hbase - 0.0005;
    n = 30;

    y = linspace(yMin, yMax, n).';
    x = xMid * ones(n, 1);

    internal = [x, y];
end


function Lconv = calcular_longitud_convectiva(nodes, elems)

    BID_CONVECCION = 2;

    Lconv = 0.0;

    for i = 1:size(elems, 1)

        if elems(i, 3) ~= BID_CONVECCION
            continue;
        end

        n1 = elems(i, 1);
        n2 = elems(i, 2);

        p1 = nodes(n1, :);
        p2 = nodes(n2, :);

        Lconv = Lconv + hypot(p2(1) - p1(1), p2(2) - p1(2));
    end
end


function datos = leer_salida_polinbe(filename)

    datos.nodes = [];
    datos.internal = [];

    fid = fopen(filename, 'r');

    if fid == -1
        error('No se pudo abrir el fichero %s', filename);
    end

    modo = "";

    while ~feof(fid)

        linea = strtrim(fgetl(fid));

        if contains(linea, '%_node_')
            modo = "node";
            continue;
        end

        if contains(linea, '%_intp_')
            modo = "intp";
            continue;
        end

        nums = sscanf(linea, '%f');

        if isempty(nums)
            continue;
        end

        if modo == "node" && numel(nums) >= 5
            datos.nodes(end + 1, :) = nums(1:5).';
        elseif modo == "intp" && numel(nums) >= 6
            datos.internal(end + 1, :) = nums(1:6).';
        end
    end

    fclose(fid);
end


function Tmedia = calcular_temperatura_media_base(datos)

    nodos = datos.nodes;

    x = nodos(:, 2);
    y = nodos(:, 3);
    T = nodos(:, 4);

    tol = 1e-10;
    idxBase = abs(y) < tol;

    xBase = x(idxBase);
    TBase = T(idxBase);

    [xBase, orden] = sort(xBase);
    TBase = TBase(orden);

    [xUnico, ~, ic] = unique(xBase);
    TUnico = accumarray(ic, TBase, [], @mean);

    if numel(xUnico) < 2
        Tmedia = mean(TBase);
    else
        Tmedia = trapz(xUnico, TUnico) / (max(xUnico) - min(xUnico));
    end
end


function analizar_y_graficar_material(nombreMaterial, Lconv, Tbase, Ta, rootDir, carpetaFiguras)

    modeloLineal = ajustar_lineal(Lconv, Tbase);
    modeloCuadratico = ajustar_cuadratico(Lconv, Tbase);
    modeloInverso = ajustar_inverso(Lconv, Tbase, Ta);
    modeloExp = ajustar_exponencial(Lconv, Tbase, Ta);

    modelos = {modeloLineal, modeloCuadratico, modeloInverso, modeloExp};

    rmse = zeros(numel(modelos), 1);

    for i = 1:numel(modelos)
        rmse(i) = modelos{i}.rmse;
    end

    [~, iBest] = min(rmse);
    mejor = modelos{iBest};

    fprintf('\nResultados regresión - %s\n', nombreMaterial);

    for i = 1:numel(modelos)
        fprintf('%s: R2 = %.5f | RMSE = %.5f ºC | coeficientes = %s\n', ...
            modelos{i}.nombre, modelos{i}.R2, modelos{i}.rmse, mat2str(modelos{i}.p, 8));
    end

    fprintf('Modelo seleccionado: %s\n', mejor.nombre);

    Lfit = linspace(min(Lconv), max(Lconv), 300).';
    Tfit = mejor.pred(Lfit);

    guardar_grafica_modelo_seleccionado(nombreMaterial, Lconv, Tbase, Lfit, Tfit, mejor, carpetaFiguras);
    guardar_grafica_residuos(nombreMaterial, Lconv, Tbase, mejor, carpetaFiguras);
    guardar_graficas_individuales_modelos(nombreMaterial, Lconv, Tbase, Lfit, modelos, carpetaFiguras);

    Tmin = min(Tbase);
    Tmax = max(Tbase);

    if Tmax < 105
        Lcrit = NaN;
        textoUmbral = "Todos los casos evaluados quedan por debajo de 105 ºC.";
        fprintf('%s\n', textoUmbral);
    elseif Tmin > 105
        Lcrit = NaN;
        textoUmbral = "Ningún caso evaluado queda por debajo de 105 ºC.";
        fprintf('%s\n', textoUmbral);
    else
        Lcrit = buscar_umbral_105(mejor, min(Lconv), max(Lconv));
        textoUmbral = sprintf("Longitud convectiva estimada para Tbase = 105 ºC: %.6f m", Lcrit);
        fprintf('%s\n', textoUmbral);
    end

    coeficientes = strings(numel(modelos), 1);

    for i = 1:numel(modelos)
        coeficientes(i) = mat2str(modelos{i}.p, 10);
    end

    seleccionado = false(numel(modelos), 1);
    seleccionado(iBest) = true;

    tablaResumen = table( ...
        string({modelos{1}.nombre; modelos{2}.nombre; modelos{3}.nombre; modelos{4}.nombre}), ...
        [modelos{1}.R2; modelos{2}.R2; modelos{3}.R2; modelos{4}.R2], ...
        [modelos{1}.rmse; modelos{2}.rmse; modelos{3}.rmse; modelos{4}.rmse], ...
        coeficientes, ...
        seleccionado, ...
        'VariableNames', {'Modelo', 'R2', 'RMSE_C', 'Coeficientes', 'Seleccionado'});

    writetable(tablaResumen, fullfile(rootDir, ['metricas_regresion_', nombreMaterial, '.csv']));

    tablaUmbral = table( ...
        string(nombreMaterial), ...
        string(mejor.nombre), ...
        Lcrit, ...
        string(textoUmbral), ...
        'VariableNames', {'Material', 'ModeloSeleccionado', 'Lcrit_m', 'Comentario'});

    writetable(tablaUmbral, fullfile(rootDir, ['umbral_105_', nombreMaterial, '.csv']));
end


function guardar_grafica_modelo_seleccionado(nombreMaterial, Lconv, Tbase, Lfit, Tfit, mejor, carpetaFiguras)

    figure('Color', 'w');
    hold on;

    scatter(Lconv, Tbase, 45, ...
        'MarkerFaceColor', [0.000 0.447 0.741], ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.7);

    plot(Lfit, Tfit, 'r-', 'LineWidth', 2.5);

    yline(105, 'k--', 'LineWidth', 2);

    aplicar_estilo_grafica();

    xlabel('Longitud convectiva [m]', 'Color', 'k', 'FontSize', 13);
    ylabel('Temperatura promedio en la base [ºC]', 'Color', 'k', 'FontSize', 13);
    title(['Regresión seleccionada - ', nombreMaterial], ...
        'Color', 'k', ...
        'FontSize', 14, ...
        'FontWeight', 'bold');

    leg = legend( ...
        'Casos BEM', ...
        ['Ajuste: ', mejor.nombre], ...
        'T máxima nominal 105 ºC', ...
        'Location', 'best');

    leg.TextColor = 'k';
    leg.Color = 'w';
    leg.EdgeColor = 'k';

    nombreFigura = fullfile(carpetaFiguras, ['regresion_seleccionada_', nombreMaterial, '.png']);
    exportgraphics(gcf, nombreFigura, 'Resolution', 300);
end


function guardar_grafica_residuos(nombreMaterial, Lconv, Tbase, mejor, carpetaFiguras)

    figure('Color', 'w');
    hold on;

    residuo = Tbase - mejor.pred(Lconv);

    scatter(Lconv, residuo, 45, ...
        'MarkerFaceColor', [0.000 0.447 0.741], ...
        'MarkerEdgeColor', 'k', ...
        'LineWidth', 0.7);

    yline(0, 'r--', 'LineWidth', 2);

    aplicar_estilo_grafica();

    xlabel('Longitud convectiva [m]', 'Color', 'k', 'FontSize', 13);
    ylabel('Residuo [ºC]', 'Color', 'k', 'FontSize', 13);
    title(['Residuos del modelo seleccionado - ', nombreMaterial], ...
        'Color', 'k', ...
        'FontSize', 14, ...
        'FontWeight', 'bold');

    nombreFigura = fullfile(carpetaFiguras, ['residuos_', nombreMaterial, '.png']);
    exportgraphics(gcf, nombreFigura, 'Resolution', 300);
end


function guardar_graficas_individuales_modelos(nombreMaterial, Lconv, Tbase, Lfit, modelos, carpetaFiguras)

    for i = 1:numel(modelos)

        modelo = modelos{i};
        Tfit = modelo.pred(Lfit);

        figure('Color', 'w');
        hold on;

        scatter(Lconv, Tbase, 45, ...
            'MarkerFaceColor', [0.000 0.447 0.741], ...
            'MarkerEdgeColor', 'k', ...
            'LineWidth', 0.7);

        plot(Lfit, Tfit, 'r-', 'LineWidth', 2.5);

        yline(105, 'k--', 'LineWidth', 2);

        aplicar_estilo_grafica();

        xlabel('Longitud convectiva [m]', 'Color', 'k', 'FontSize', 13);
        ylabel('Temperatura promedio en la base [ºC]', 'Color', 'k', 'FontSize', 13);

        titulo = sprintf('%s - %s | R2 = %.4f | RMSE = %.4f ºC', ...
            nombreMaterial, modelo.nombre, modelo.R2, modelo.rmse);

        title(titulo, ...
            'Color', 'k', ...
            'FontSize', 13, ...
            'FontWeight', 'bold');

        leg = legend( ...
            'Casos BEM', ...
            ['Ajuste: ', modelo.nombre], ...
            'T máxima nominal 105 ºC', ...
            'Location', 'best');

        leg.TextColor = 'k';
        leg.Color = 'w';
        leg.EdgeColor = 'k';

        nombreModelo = matlab.lang.makeValidName(modelo.nombre);
        nombreFigura = fullfile(carpetaFiguras, ...
            ['regresion_', nombreMaterial, '_', nombreModelo, '.png']);

        exportgraphics(gcf, nombreFigura, 'Resolution', 300);
    end
end


function aplicar_estilo_grafica()

    grid on;
    box on;

    ax = gca;
    ax.Color = 'w';
    ax.XColor = 'k';
    ax.YColor = 'k';
    ax.GridColor = [0.75 0.75 0.75];
    ax.GridAlpha = 0.45;
    ax.LineWidth = 1.2;
    ax.FontSize = 12;
end


function modelo = ajustar_lineal(L, T)

    p = polyfit(L, T, 1);

    pred = @(x) polyval(p, x);

    [R2, rmse] = metricas(T, pred(L));

    modelo.nombre = 'Lineal';
    modelo.p = p;
    modelo.pred = pred;
    modelo.R2 = R2;
    modelo.rmse = rmse;
end


function modelo = ajustar_cuadratico(L, T)

    p = polyfit(L, T, 2);

    pred = @(x) polyval(p, x);

    [R2, rmse] = metricas(T, pred(L));

    modelo.nombre = 'Cuadrático';
    modelo.p = p;
    modelo.pred = pred;
    modelo.R2 = R2;
    modelo.rmse = rmse;
end


function modelo = ajustar_inverso(L, T, Ta)

    fun = @(p, x) Ta + p(1) ./ (x + p(2));

    obj = @(p) mean((fun(p, L) - T).^2);

    p0 = [mean((T - Ta) .* L), 0.01];

    p = fminsearch(obj, p0, optimset('Display', 'off'));

    pred = @(x) fun(p, x);

    [R2, rmse] = metricas(T, pred(L));

    modelo.nombre = 'Inverso asintótico';
    modelo.p = p;
    modelo.pred = pred;
    modelo.R2 = R2;
    modelo.rmse = rmse;
end


function modelo = ajustar_exponencial(L, T, Ta)

    Y = max(T - Ta, 1e-6);

    pLog = polyfit(L, log(Y), 1);

    a = exp(pLog(2));
    b = pLog(1);

    pred = @(x) Ta + a * exp(b * x);

    [R2, rmse] = metricas(T, pred(L));

    modelo.nombre = 'Exponencial asintótico';
    modelo.p = [a, b];
    modelo.pred = pred;
    modelo.R2 = R2;
    modelo.rmse = rmse;
end


function [R2, rmse] = metricas(Treal, Tpred)

    err = Tpred - Treal;

    ssRes = sum(err.^2);
    ssTot = sum((Treal - mean(Treal)).^2);

    R2 = 1 - ssRes / ssTot;
    rmse = sqrt(mean(err.^2));
end


function Lcrit = buscar_umbral_105(modelo, Lmin, Lmax)

    f = @(L) modelo.pred(L) - 105;

    if f(Lmin) * f(Lmax) > 0
        Lcrit = NaN;
        return;
    end

    Lcrit = fzero(f, [Lmin, Lmax]);
end