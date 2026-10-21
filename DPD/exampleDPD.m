
%% CONFIGURACIÓN

projectRoot = fileparts(fileparts(mfilename('fullpath')));
addpath(fullfile(projectRoot, 'DPD'));

medidas = fullfile(projectRoot, 'measurements', ...
    'experiment20260429T134032_xy.mat');

entradaOnline = medidas;
P = 310;

%% OFFLINE + ONLINE

% Entrenar los tres modelos desde cero.
modelFile = trainDPD(medidas, P);

% Si el OFFLINE ya está hecho, comentar la línea anterior
% y descomentar esta para reutilizar los modelos guardados.
% modelFile = fullfile(projectRoot, 'results', 'dpd_sweep', ...
%     'tp98e6f7ac_f3cd_44c7_9ce3_fd8619319748', 'models_P0310.mat');

% Generar las tres señales DPD.
[outputFile, dpd] = runDPDOnline(modelFile, entradaOnline);

%% SEÑALES DPD

xGMP  = dpd(1).yvalmod;
xPNIQ = dpd(2).yvalmod;
xPNNN = dpd(3).yvalmod;

%% RESULTADO

fprintf('Modelos: %s\n', modelFile);
fprintf('Señales guardadas: %s\n', outputFile);
