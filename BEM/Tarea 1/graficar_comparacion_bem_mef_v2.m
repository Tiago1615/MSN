clear; clc; close all;

fMEF = 'Puntos_internos_MEF.txt';

fConOut = 'poconbe.txt';
fLinOut = 'polinbe.txt';

% Leer MEF y eliminar puntos de contorno
mef = readmatrix(fMEF, 'FileType', 'text', 'CommentStyle', '%');
mef = mef(~any(isnan(mef), 2), :);

yMEF_all = mef(:, 1);
TMEF_all = mef(:, 2);

idx_int = yMEF_all > min(yMEF_all) + 1e-12 & ...
          yMEF_all < max(yMEF_all) - 1e-12;

yMEF = yMEF_all(idx_int);
TMEF = TMEF_all(idx_int);

% Leer resultados BEM
[yCon, TCon] = read_bem_internal(fConOut);
[yLin, TLin] = read_bem_internal(fLinOut);

% Gráfica temperaturas
figure('Color', 'w');
hold on;

plot(TMEF, yMEF, 'k-', 'LineWidth', 2);

legend_entries = {'MEF / COMSOL'};

if ~isempty(yCon)
    plot(TCon, yCon, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4);
    legend_entries{end + 1} = 'PoConBE';
end

if ~isempty(yLin)
    plot(TLin, yLin, 's-', 'LineWidth', 1.2, 'MarkerSize', 4);
    legend_entries{end + 1} = 'PoLinBE';
end

grid on; box on;
xlabel('Temperatura [ºC]');
ylabel('Altura en el plano medio de la aleta central [m]');
title('Temperatura en puntos internos: MEF vs PoConBE vs PoLinBE');
legend(legend_entries, 'Location', 'best');

% Gráfica errores
figure('Color', 'w');
hold on;

legend_err = {};

if ~isempty(yCon)
    TCon_i = interp1(yCon, TCon, yMEF, 'linear', 'extrap');
    errCon = TCon_i - TMEF;

    plot(errCon, yMEF, 'o-', 'LineWidth', 1.2, 'MarkerSize', 4);
    legend_err{end + 1} = 'PoConBE';

    fprintf('PoConBE: max abs = %.6g ºC | RMS = %.6g ºC\n', ...
        max(abs(errCon)), sqrt(mean(errCon.^2)));
end

if ~isempty(yLin)
    TLin_i = interp1(yLin, TLin, yMEF, 'linear', 'extrap');
    errLin = TLin_i - TMEF;

    plot(errLin, yMEF, 's-', 'LineWidth', 1.2, 'MarkerSize', 4);
    legend_err{end + 1} = 'PoLinBE';

    fprintf('PoLinBE: max abs = %.6g ºC | RMS = %.6g ºC\n', ...
        max(abs(errLin)), sqrt(mean(errLin.^2)));
end

xline(0, 'k--');

grid on; box on;
xlabel('T_{BEM} - T_{MEF} [ºC]');
ylabel('Altura [m]');
title('Error frente a MEF');

if ~isempty(legend_err)
    legend(legend_err, 'Location', 'best');
end

% --------------------------------------------------------
function [y,T] = read_bem_internal(filename)

y = [];
T = [];

if ~isfile(filename)
    warning('No existe %s', filename);
    return;
end

fid = fopen(filename,'r');

if fid==-1
    warning('No se pudo abrir %s', filename);
    return;
end

modo = 0;
datos = [];

while ~feof(fid)

    linea = strtrim(fgetl(fid));

    if contains(linea,'%_intp_')
        modo = 1;
        continue;
    end

    if modo==1

        nums = sscanf(linea,'%f');

        if numel(nums) >= 4
            % formato:
            % id x y T qx qy
            datos(end+1,:) = nums(:).';
        end
    end
end

fclose(fid);

if isempty(datos)
    warning('%s no contiene puntos internos.',filename);
    return;
end

y = datos(:,3);
T = datos(:,4);

[y,idx] = sort(y);
T = T(idx);

[y,ia] = unique(y,'stable');
T = T(ia);

end