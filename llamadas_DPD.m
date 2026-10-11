
%% CONFIGURACIÓN

medidasEntrenamiento = 'measurements/experiment20260429T134032_xy.mat';
entradaOnline = 'measurements/experiment20260429T134032_xy.mat';

P = 340;

entrenar = true;  % true:  OFFLINE + ONLINE
                  % false: SOLO ONLINE

modeloGuardado = ['results\dpd_offline\' ...
    'P340_20261011_210433_267_ed6696d7-1826-412c-b689-f7f61a53abfa\' ...
    'models.mat'];

%% OFFLINE

if entrenar
    modelFile = trainDPDOffline(medidasEntrenamiento, P);
else
    modelFile = modeloGuardado;
end

%% ONLINE

outputFile = runDPDOnline(modelFile, entradaOnline);

%% EXTRAER SEÑALES

S = load(outputFile, 'dpd');

% comprobar amplitud, CFR, etc ???
xGMP  = S.dpd(1).yvalmod;
xPNIQ = S.dpd(2).yvalmod;
xPNNN = S.dpd(3).yvalmod;

%% RESUMEN

fprintf('\nModelo utilizado:\n%s\n', modelFile);
fprintf('\nSeñales guardadas en:\n%s\n', outputFile);

fprintf('\nMuestras generadas:\n');
fprintf('GMP-DOMP: %d\n', numel(xGMP));
fprintf('PN-IQ-GMP: %d\n', numel(xPNIQ));
fprintf('PNNN: %d\n', numel(xPNNN));
