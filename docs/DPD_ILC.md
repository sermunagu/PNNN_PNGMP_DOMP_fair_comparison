# DPD offline por imitación del ILC

> Organización actual (14/10/2026): ejecutores PA en `PA/` y DPD en `DPD/`;
> `runPNIQAblation` y `runPNNNPhaseExperiment` están en `DPD/experiments/`.
> Antes de usar ejemplos: `addpath('PA','DPD','PA/experiments','DPD/experiments')`.
> Las narraciones históricas conservan sus nombres; el [mapeo actual](GUIA_REPOSITORIO_DPD.md#entradas-reorganizadas-14102026)
> indica sus sustitutos, sin cambiar medidas, modelos ni resultados.

> Estado actualizado a 13/10/2026: el barrido offline real P20:10:500 está
> terminado y selecciona P310. También se completaron la ablación
> P100/P310/P500 y el experimento PNNN de cero estructural P310.
> Los apartados anteriores que los presentan como pendientes describen
> etapas históricas. La validación RF sigue pendiente; véase la
> [guía actual](GUIA_REPOSITORIO_DPD.md).

Los tests de desarrollo mencionados en este documento se conservan fuera
del árbol principal; la guía indica su copia recuperable y cómo restaurarlos.

## OFFLINE / ONLINE

Desde la raíz del repositorio, tras `addpath('PA','DPD','PA/experiments','DPD/experiments')`:

```matlab
modelFile = trainDPD('measurements/nueva_medida.mat', 340);
outputFile = runDPDOnline(modelFile, 'measurements/entrada_online.mat');
```

OFFLINE carga directamente `x` e `y`. Cada llamada ajusta de nuevo los tres
modelos, conservando el centrado DC, la identificación, DOMP, N12/sigmoide,
selección de épocas, poda y fine-tuning existentes. No reutiliza entrenamientos.

Una carpeta nueva creada con `tempname` bajo `results/dpd_offline/` contiene:

- `models.mat`: los tres modelos en `bundle`, configuración, split, métricas,
  archivo de origen y `fs` si existe.
- `predictions.mat`: predicciones directas del ajuste.

ONLINE carga el `bundle`, lee `x` y llama a `applyDPDModels`, sin entrenar.
Guarda un MAT nuevo junto al modelo, con `dpd(k).yvalmod` y `dpd(k).modeltype`
en orden GMP-DOMP, PN-IQ-GMP, PNNN. Tras reiniciar MATLAB, usar la ruta del
`models.mat` guardado. Los paquetes existentes siguen siendo compatibles.

Las rutas relativas las resuelve MATLAB. Se usa el formato conocido del tutor;
archivo o campos incorrectos producen los errores normales de MATLAB.
No hay helpers de carga, hashes de medidas, UUID ni comprobaciones de calibración.
El investigador debe usar señales compatibles en frecuencia, escala y periodicidad.
No se escala, recorta ni restaura DC; el NMSE offline no es una medida RF.

Se conservan los tests sintéticos ligeros de entrenamiento y recarga. El tercer
argumento opcional `cfg` de OFFLINE permite reducir épocas y usar una carpeta
temporal en esos tests; las dos llamadas normales mantienen el protocolo original.

Tras la simplificación pasan los tests sintéticos y la recarga en otra sesión
MATLAB. El P340 existente reproduce exactamente, señal por señal y con el mismo
tipo numérico, su archivo ONLINE anterior (error relativo 0 en los tres modelos).
No se repitió su entrenamiento. `checkcode` no informa diagnósticos.

## Barrido de modelado PA y reanudación

El barrido canónico sigue usando por defecto la captura `_forward_xy.mat`;
no es el barrido DPD ni cambia el significado físico de x/y.

```matlab
sweep = runPASweep(); % nueva carpeta, no reutiliza otros experimentos
sweep = runPASweep('results/parameter_sweep/carpeta_de_esta_ejecucion');
comparison = comparePA(sweep.selection.selectedParameters, sweep);
```

La segunda llamada es una alternativa para reanudar, no otro entrenamiento.
La carpeta se imprime antes de entrenar. `experiment.mat` guarda configuración,
señales preprocesadas y split; los MAT lineal, Ridge, denso y
`pnnn_target_XXXX.mat` guardan cada etapa. Solo se calculan etapas ausentes.
Cada PNNN parte de la misma fuente densa, no de la poda anterior.
No cambie los archivos de esa carpeta para reutilizarla con otra captura.

Los nuevos checkpoints usan variables directas (`linear`, `fixedResults`,
`denseSource`, `point`, `result`), sin firmas ni versiones de protocolo.
No se migran ni sobrescriben checkpoints históricos: para ejecutarlos con el
formato anterior debe usarse la revisión histórica del código.
Para regenerar figuras de una ejecución nueva se carga
`complexity_sweep.mat` y se pasa `S.result` a `comparePA`.
Reanudar vuelve a exportar las tres figuras del barrido, sin repetir ajustes.

Las métricas y checkpoints se guardan antes de exportar. No hay preflight:
`matlab2tikz` solo se necesita al producir figuras FIG/PNG/TikZ. Un error gráfico
no borra el entrenamiento; se puede reanudar la misma carpeta.
El tercer argumento opcional `cfg` permite pruebas sintéticas reducidas,
como en OFFLINE; el uso normal mantiene el preset científico.

## Barrido DPD de complejidad (imitación del ILC)

```matlab
sweep = runDPDSweep('measurements/experiment20260429T134032_xy.mat', 20:10:500);
% Si se interrumpe, usar la carpeta impresa al comenzar:
sweep = runDPDSweep('results/dpd_sweep/carpeta_de_esta_ejecucion');
% Cualquier P calculado se puede usar en otra sesión, sin entrenamiento:
outputFile = runDPDOnline(fullfile(sweep.resultDirectory, 'models_P0200.mat'), 'measurements/entrada_online.mat');
```

La primera llamada crea una ejecución nueva. La segunda reanuda explícitamente
esa ejecución; se distingue una carpeta real de un archivo mediante `isfolder`.
La captura es obligatoria: no se hereda la captura PA por defecto. En DPD
`x=u`, `y=x_ILC`, no la salida física del amplificador. El preprocesado y ajuste
son los de `trainDPD`, sin modificar los scripts ni resultados PA.

### Qué se reutiliza del barrido PA

Se siguieron las dependencias de `main_sweep_and_comparison`,
`run_parameter_sweep` y `run_selected_comparison`, incluidos población,
diccionarios, DOMP, LS/Ridge, dataset neuronal, selección de épocas, poda,
normalización, predicción, NMSE, FLOPs y exportación. El nuevo coordinador
reutiliza las funciones científicas; no copia sus algoritmos.

- **Población GMP:** `regressorManager.buildPopulation` fija orden 13, lags 0:10 y separación de
  envolvente 1:2. Son 143 términos principales + 264 retardados + 264 adelantados.
  Se añaden x, conj(x) y abs(x); al eliminar el x duplicado quedan **673**.
  `buildGMPManager` conserva la canonicalización y el orden.
- **Dos rutas diferentes:** `fitLinearModels` comparte una población entre
  `fitComplexGMP` y `fitPNIQGMP`. Cada uno llama una vez a
  `selectDOMPSupport`, sobre identificación, hasta max(P)/2. PN-IQ usa su
  diccionario real de 1285 características (2·673 menos 61 columnas Q nulas),
  objetivo rotado y su propia selección; no hereda el soporte GMP.
- **Cada prefijo:** para P se reajustan P/2 columnas mediante `lsqminnorm`.
  La normalización vigente es **pico por columna y pico de salida**, no norma 2
  ni RMS. GMP ajusta coeficientes complejos; PN-IQ ajusta I y Q separadamente.
  Los coeficientes físicos desnormalizados generan predicciones, conservando
  la rotación/restauración PN-IQ. La métrica de rango usa los coeficientes
  normalizados (módulo complejo GMP, máximo absoluto conjunto I/Q PN-IQ).
- **Una fuente PNNN:** `trainDensePNNN` selecciona épocas con
  train/validación internos de identificación y duración de fine-tuning con
  P200. Reinicializa desde la misma semilla y entrena definitivamente con toda
  identificación, recomputando allí el z-score histórico. Son **dos llamadas
  trainnet** (provisional y definitiva), no dos fuentes para podar ni una por P.
  Se conservan 84 entradas, fc1(12), sigmoide, fcOut(2), semillas y normalización.
- **Poda independiente:** cada `fitSparsePNNN` parte de esa misma
  densa definitiva, con poda global por magnitud, 14 biases protegidos y
  fine-tuning de máscara fija (incluida máscara de gradientes). Nunca se encadena
  P20→P30. Pesos permitidos por máscara + biases = P.
- **Costes:** se reutilizan `countModelOperations`, `countModelFLOPs`,
  `countSparsePNNNFLOPs` y `summarizeTrainableParameters`. Adición/multiplicación
  real valen 1/1, compleja 2/6; operaciones especiales se cuentan aparte.
  PNNN reporta FLOPs **ideales sparse**, no tiempo MATLAB ni matriz comprimida.
- **Ridge suplementario:** `fitFixedRidge` reajusta los soportes ya
  seleccionados con los lambdas fijos existentes; no añade DOMP. Se añade su
  NMSE de evaluación. Los bundles principales siguen teniendo tres modelos.
- **Persistencia y gráficos:** se conserva el planteamiento de checkpoints
  explícitos del PA. El PA guarda métricas/predicciones lineales, densa, puntos,
  selección, CSV/MAT y tres curvas globales; el punto seleccionado añade espectros
  Welch compartidos, tiempo, AM/PM y error de fase. El nuevo DPD guarda además
  modelos compactos por P y evaluación complementaria; no llama al comparador
  PA ni modifica Welch. Reutiliza `plotSweepPaperFigure`, `getIEEEPaperStyle` y
  `exportPaperFigure`, sin el límite PA de −30 dB ni la etiqueta `Validation NMSE`.

### Artefactos y recuperación

Cada ejecución usa una carpeta nueva bajo `results/dpd_sweep/`:

| Archivo | Contenido |
|---|---|
| `experiment.mat` | `cfg`, x/y preprocesados, `split` con índices de evaluación, `provenance` (archivo, medias, semántica, fs si existe). No necesita releer medidas al reanudar. |
| `linear_sweep.mat` | `linear` (tablas, rutas y predicciones por P), `linearModels` (coeficientes y descriptores máximos compartidos). |
| `fixed_lambda_linear_sweep.mat/.csv` | `fixedResults`, suplementarios con NMSE identificación/evaluación/full, lambda y costes. |
| `sweep_dense_source.mat` | `denseSource` definitivo, normalización y épocas seleccionadas. |
| `pnnn_target_XXXX.mat` | `point`: modelo, máscaras, fila de métricas y predicción completa. |
| `models_PXXXX.mat` | `bundle` autocontenido compatible con ONLINE, exactamente tres modelos, config, split, métricas y procedencia. |
| `dpd_complexity_sweep.mat/.csv` | `sweep` y su tabla `results`; el MAT incluye configuración, selección/razón de no selección, metadatos, rutas y tabla Ridge. |
| `operating_point_selection*.csv` | Diagnóstico y sensibilidad, solo si se encuentra un punto admisible. |
| `rf_measurements.csv` | Model, P, ModelFile, CaptureFile, RFMeasurementFile vacío para asociar futuras medidas reales; se conserva al reanudar. |

Las columnas principales son `Model`, `TargetRealParameters`,
`ActualRealParameters`, `SelectedLambda`, `InternalValidationNMSEdB`
(NaN cuando no aplica), `IdentificationNMSEdB`, `FullSignalNMSEdB`,
`FLOPsPerSample`, `ActiveWeights`, `ActiveBiases`, `WeightSparsityPercent`,
`FineTuneEpochs`, `MaxAbsRealParameter` y `EvaluationNMSEdB`.

Cada GMP ONLINE almacena solo P/2 regresores y coeficientes. PN-IQ guarda el
prefijo de `featureMap` y únicamente sus regresores fuente únicos, descriptores
correspondientes y `selectedColumns` recalculados con el **nuevo desplazamiento
Q local**. No se evalúa el diccionario de P máximo con ceros sobrantes.
`applyDPDModels` y `runDPDOnline` no han cambiado.

Se exportan nueve figuras FIG/PNG/TikZ: NMSE identificación, evaluación y full
frente a P y FLOPs; FLOPs frente a P; sparsity PNNN frente a P; rango máximo
frente a P. Las curvas frente a P incluyen Ridge cuando procede. El rango
lineal y el neuronal usan parametrizaciones diferentes; no son ganancias físicas
directamente comparables. Las métricas se guardan antes de exportar figuras.
Para repetir solo exportación, cargar el MAT y ejecutar
`plotDPDSweepFigures(S.sweep)` con `toolbox/` en el path; no se cargan modelos.
En las leyendas se abrevia Complex GMP-DOMP como GMP-DOMP para dejar espacio
a las referencias Ridge; las tablas y modelos conservan sus nombres completos.

### Evaluación y selección: límites científicos

`EvaluationNMSEdB` usa el complemento exacto de identificación y la función común
`nmseComplexDb`. No se utiliza para ajustar pesos/soportes ni seleccionar épocas,
semillas o hiperparámetros. Sin embargo, **no es un test externo independiente**:
`selectPeakSegment` sitúa el segmento de identificación atendiendo al máximo de y de
la captura completa; el centrado DC usa la captura completa; los retardos
periódicos pueden atravesar los bordes entre subconjuntos. No se cambió nada de
ello. La validación interna PNNN está contenida en identificación.

Se aplica `selectOperatingPoint` a las curvas de evaluación con ventana y
tolerancia originales. No se presupone P340. Si falta una ventana completa,
estabilización conjunta o monotonía de FLOPs, se guardan todas las curvas con
`selection=[]` y `selectionReason`; no se fuerza un P. Usar evaluación para
elegir P la convierte en conjunto de desarrollo/selección, no confirmación final.
Una confirmación futura requiere otra captura independiente, o una partición
predefinida sin selección por y y con guardas temporales justificadas por memoria,
preprocesado fijado con datos de entrenamiento y tratamiento explícito de bordes.

Para un artículo RF faltan medidas del PA con las tres señales DPD y referencias
comparables, sincronización, potencia/ganancia/CFR documentados y ACLR/ACPR/EVM
cuando procedan. No se calculan métricas RF a partir de imitación ILC ni se
atribuye mejora de linealización a estas curvas offline.

### Pruebas pequeñas

`run_dpd_parameter_sweep_test(folder)` usa dos capturas sintéticas de 512 muestras,
P20/30/40, arquitectura original y dos épocas históricas. Instrumenta DOMP y
trainnet, verifica poda desde la misma densa, parámetros/FLOPs/NMSE, compacción
incluido Q, OFFLINE independiente, ONLINE y reanudación total/parcial. Conserva
un fixture para `run_dpd_online_fresh_session_test(folder)` en otra sesión MATLAB.
El tercer argumento cfg permite `cfg.dpd.sweepResultsRoot` y
`cfg.dpd.exportFigures=false` para estas pruebas; no altera el uso normal.
Las tolerancias son relativas L2 1e-10 para doble y 1e-6 para predicción neuronal
single, comparadas por familia para evitar que concatenar PNNN convierta los
lineales a single. Máscaras y soportes se comparan exactamente.

Validación del 12/10/2026 en MATLAB R2026a: test de barrido aprobado (~45 s),
recarga en sesión nueva aprobada, nueve exportaciones FIG/PNG/TikZ y apertura
visible verificadas (~50 s). También pasan los tests existentes
`run_linear_complexity_sweep_test`, `run_fixed_lambda_linear_sweep_test`,
`run_pnnn_shared_dense_sweep_test` y `run_operating_point_selection_test`.
Los cuatro archivos MATLAB nuevos pasan `checkcode` sin diagnósticos.
La prueba gráfica `run_dpd_sweep_presentation_test(folder)` verifica datos,
CSV y ausencia de ajustes durante exportación. No se ejecutó el barrido real,
no se reentrenó P340 y no hubo transmisión RF.

## Experimentos históricos

La recuperación P340 del 9 de octubre reajustó los modelos lineales sobre
soportes guardados y repitió únicamente poda/fine-tuning desde la densa.
No se repitieron sus 1228 épocas. Los archivos históricos no se migran.
`testDPD.m` ya estaba eliminado del working tree al iniciar esta limpieza;
no es una llamada operativa vigente. Para inferir use el `models.mat` existente
con ONLINE. La captura DPD es `experiment20260429T134032_xy.mat`:
`x=u`, `y=x_ILC`; no debe confundirse con el modelado PA.

El registro de entrada es periódico. ONLINE centra su propia media y conserva
la normalización entrenada. No restaura DC ni adapta pico, CFR o límites DAC.
La compatibilidad física de las nuevas señales debe comprobarla el investigador.

## Auditoría de P=340

| Familia | Grados reales activos del modelo | Almacenamiento numérico principal |
|---|---|---|
| Complex GMP | 170 coeficientes complejos = 340 reales | 170 complejos físicos; además coeficientes de comparación y metadatos |
| PN-IQ | 170 características, 170 coeficientes I + 170 Q | 340 reales físicos; además comparación y metadatos |
| PNNN | 326 pesos en máscara + 14 biases protegidos | 1046 posiciones densas (1032 pesos + 14 biases), incluidas las podadas a cero |

PNNN usa 84 entradas, fc1(12), sigmoide y fcOut(2); la sigmoide no tiene
parámetros entrenables. Activo significa grado permitido por soporte/máscara,
no `nnz` del valor numérico: un bias activo puede valer cero. No se cuentan
estadísticas mu/sigma, regresores, índices, máscaras ni constantes de activación
como parámetros ajustados. El criterio compara grados ajustables, no bytes de
almacenamiento ni capacidad funcional idéntica. Las columnas redundantes o
estructuralmente constantes pueden reducir grados efectivos: no se altera aquí
el criterio histórico ni se elimina ninguna característica.

FLOPs es otra magnitud: el checkpoint reporta 1817/1107/840 por muestra,
respectivamente. La cifra sparse de PNNN supone aprovechar pesos podados;
MATLAB almacena matrices densas y no garantiza ese coste real de ejecución.

Limitación numérica detectada en la prueba sintética: el canal imaginario del
tap actual, teóricamente cero tras normalizar fase, puede tener sigmaX~1e-17.
La regla histórica solo sustituye sigma exactamente cero; puede amplificar
redondeo y romper la equivariancia de fase numérica (1.03% L2 en este fixture).
La inferencia conserva esa regla y reproduce la evaluación original. Una
corrección futura sería fijar explícitamente ese canal a cero o umbralizar su
sigma, pero requiere aprobación y nueva validación/entrenamiento; no se aplica
silenciosamente ni se atribuye ese porcentaje al experimento real.
En el checkpoint denso real se confirma sigmaX=8.019388e-19 en la característica
15; no se ha medido aquí el efecto de rotar la captura real.

## Validación y límites

```matlab
addpath('tests');
run_dpd_offline_online_test;
run_dpd_model_reuse_test;
run_pnnn_shared_dense_sweep_test;
run_sweep_resume_test;
```

Los tests usan capturas sintéticas y pocas épocas. Comprueban soportes,
presupuestos, máscaras, guardado/recarga, ONLINE sin entrenamiento y
reanudación completa o con un único punto ausente. El test de regresión de
simplificación permite registrar una referencia antes de editar y compararla
después; no genera su propio resultado esperado durante la comprobación.

Siguen pendientes la evaluación sobre muestras externas, la decisión sobre
sigma casi nula, y la transmisión y medida RF autorizadas. El NMSE offline
sobre señal completa incluye muestras de identificación; no es un holdout
independiente ni una medida de linealización del PA.
