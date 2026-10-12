
%% CONFIGURACIÓN

projectRoot = fileparts(mfilename('fullpath'));
addpath(projectRoot);

% DPD por imitación ILC: x = entrada deseada, y = objetivo ILC.
medidasEntrenamiento = fullfile(projectRoot, 'measurements', 'experiment20260429T134032_xy.mat');
entradaOnline = medidasEntrenamiento; % Puede ser otro MAT; solo necesita x.

P = 340;

entrenar = false; % true: entrenar desde cero y guardar en una carpeta nueva.
                  % false: solo inferencia con modeloGuardado, sin entrenar.

modeloGuardado = fullfile(projectRoot, 'results', 'dpd_offline', ...
    'P340_20261011_210433_267_ed6696d7-1826-412c-b689-f7f61a53abfa', ...
    'models.mat'); % En ONLINE se usa el presupuesto del modelo guardado.

%% OFFLINE

if entrenar
    modelFile = trainDPDOffline(medidasEntrenamiento, P);
else
    modelFile = modeloGuardado;
end

%% ONLINE

[outputFile, dpd] = runDPDOnline(modelFile, entradaOnline);

%% EXTRAER SEÑALES

% Sin transmisión RF ni ajuste de amplitud, DC o CFR.
xGMP  = dpd(1).yvalmod;
xPNIQ = dpd(2).yvalmod;
xPNNN = dpd(3).yvalmod;

%% RESUMEN

fprintf('\nModelo utilizado:\n%s\n', modelFile);
fprintf('\nSeñales guardadas en:\n%s\n', outputFile);

fprintf('\nMuestras generadas:\n');
fprintf('GMP-DOMP: %d\n', numel(xGMP));
fprintf('PN-IQ-GMP: %d\n', numel(xPNIQ));
fprintf('PNNN: %d\n', numel(xPNNN));
