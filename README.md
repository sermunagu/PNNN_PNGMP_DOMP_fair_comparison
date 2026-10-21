# PA Modeling / Digital Predistortion

MATLAB: Complex GMP-DOMP, PN-IQ-GMP y PNNN. Los ejecutores de modelado
del amplificador (`PA/`) y predistorsión (`DPD/`) comparten los algoritmos
de `toolbox/`. Medidas y resultados históricos mantienen sus ubicaciones.

Desde la raíz del repositorio, con MATLAB R2026a, Deep Learning Toolbox y
Signal Processing Toolbox:

```matlab
addpath('PA','DPD','PA/experiments','DPD/experiments');
```

Las entradas añaden `config/` y las partes necesarias de `toolbox/` usando
su propia ubicación. No usar `addpath(genpath(pwd))`: incluiría la referencia
independiente `pn-iq-gmp/` con funciones homónimas. La exportación utiliza
`third_party/matlab2tikz` y produce FIG, PNG y TikZ; el PDF final se compila en Overleaf.

## PA Modeling

Aprende entrada del PA → salida medida. La captura por defecto es
`measurements/experiment20260429T134032_forward_xy.mat`.

```matlab
% NUEVO: costoso; barrido, selección y comparación del punto seleccionado.
comparison = runPA();
% Alternativa: solo un barrido NUEVO, sin la comparación posterior.
pa = runPASweep(20:10:500);
% REANUDAR una ejecución del formato actual, reutilizando checkpoints.
pa = runPASweep('results/parameter_sweep/carpeta_de_la_ejecucion');
% COMPARAR: predicciones guardadas, refit fixed Ridge y exportación de figuras.
comparison = comparePA(pa.selection.selectedParameters, pa);
% CONSULTAR: únicamente leer, sin ajustar ni exportar.
S = load('results/parameter_sweep/carpeta_de_la_ejecucion/complexity_sweep.mat');
S.result.results
% Contar candidatos, sin cargar medidas ni entrenar.
countGMPParams
```

Son alternativas, no una secuencia para ejecutar dos barridos. Atención:
`comparePA(P)` **sin segundo argumento inicia un barrido completo**; se
conserva ese contrato anterior. Pasar siempre `pa` para reutilizarlo.
Reanudar regenera tablas/figuras, por lo que tampoco es una consulta de solo lectura.
Los checkpoints PA legacy sin `experiment.mat` se conservan como referencias
de lectura, no son compatibles con la reanudación/comparación del formato actual.
El selector del histórico PA sigue dando **P340**.

`probeWideGMP(referenceCSV)`, en `PA/experiments/`, es un experimento PA
costoso P340 con selección DOMP/Ridge propia; no es la ablación DPD.

## Digital Predistortion

Aprende por imitación del ILC: **x = u**, entrada deseada; **y = x_ILC**,
predistorsión objetivo. Captura histórica: `measurements/experiment20260429T134032_xy.mat`.

```matlab
% OFFLINE NUEVO: entrena los tres modelos desde cero en una carpeta nueva.
modelFile = trainDPD('measurements/nueva_medida.mat',310);
% ONLINE: solo necesita x; no entrena ni transmite RF.
[outputFile,dpd] = runDPDOnline(modelFile,'measurements/entrada_online.mat');
% BARRIDO NUEVO: costoso; una fuente densa compartida y podas independientes.
sweep = runDPDSweep('measurements/nueva_medida.mat',20:10:500);
% REANUDAR explícitamente la misma ejecución.
sweep = runDPDSweep('results/dpd_sweep/carpeta_de_la_ejecucion');
% CONSULTAR sin entrenamiento ni exportación.
S = load('results/dpd_sweep/carpeta_de_la_ejecucion/dpd_complexity_sweep.mat');
S.sweep.results
```

`dpd(1:3).yvalmod` conserva el orden GMP-DOMP, PN-IQ-GMP, PNNN, con sus
etiquetas `modeltype`. `exampleDPD` inicia OFFLINE y ONLINE con P310; para reutilizar el modelo
guardado, comentar `trainDPD` y activar la ruta indicada en el script. El barrido DPD histórico selecciona **P310**.
Una reanudación completa no construye características PNNN ni carga la fuente densa.
Las rutas relativas suministradas se interpretan desde el directorio de trabajo;
desde otra carpeta, añadir las entradas con rutas absolutas y usar argumentos absolutos.

### Experimentos DPD

```matlab
% Ablación lineal: cinco modelos/controles; no entrena PNNN. También P100/P500.
ablation = runPNIQAblation(sweep.resultDirectory,310);
% SOLO si se desea entrenar una nueva PNNN experimental: ejecución costosa.
[correctedFile,comparison] = runPNNNPhaseExperiment(sweep.resultDirectory,310);
```

Los diagnósticos `diagnosePNIQAblation`, `diagnosePNNNPhaseNormalization` y
`auditPNNNPhaseSource` permanecen en `toolbox/` (para invocarlos directamente:
`addpath(genpath('toolbox'))`). Se conservan SVD, sensibilidad, soportes comunes
y los resultados P100/P310/P500. La variante de cero estructural es explícita;
no sustituye la normalización ni el modelo PNNN histórico.

Estos resultados son **offline**, no una validación RF. Se mantienen DC,
memoria periódica y normalizaciones; inferir requiere unidades y condiciones
compatibles con el entrenamiento. Falta validación independiente en laboratorio.

La [guía](docs/GUIA_REPOSITORIO_DPD.md) explica algoritmos, bundles, el mapeo
de nombres antiguos y la recuperación de los tests internos externos.
Los informes [PN-IQ](docs/PNIQ_GMP_ABLATION_2026-10-12.md) y
[PNNN](docs/PNNN_PHASE_NORMALIZATION_DIAGNOSIS_2026-10-12.md) conservan la evidencia.
