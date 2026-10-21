# Guía del repositorio DPD

Estado contrastado con código y archivos locales el **14/10/2026**. Raíz: `C:\Sergi\Investigacion\Códigos\PNNN_PNGMP_DOMP_fair_comparison`. Hay WIP sin commit: GitHub no es una copia completa del estado descrito aquí.

### Entradas reorganizadas (14/10/2026)

Desde la raíz: `addpath('PA','DPD','PA/experiments','DPD/experiments')`.
No hay wrappers con los nombres antiguos; los campos internos de modelos y checkpoints no cambian.

| Nombre anterior en la raíz | Ubicación actual |
|---|---|
| `main_sweep_and_comparison.m` | `PA/runPA.m` |
| `run_parameter_sweep.m` | `PA/runPASweep.m` |
| `run_selected_comparison.m` | `PA/comparePA.m` |
| `count_dense_parameters.m` | `PA/countGMPParams.m` |
| `run_widely_linear_gmp_probe.m` | `PA/experiments/probeWideGMP.m` |
| `trainDPDOffline.m` | `DPD/trainDPD.m` |
| `runDPDParameterSweep.m` | `DPD/runDPDSweep.m` |
| `runDPDOnline.m` | `DPD/runDPDOnline.m` |
| `llamadas_DPD.m` | `DPD/exampleDPD.m` |
| `runPNIQAblation.m` | `DPD/experiments/runPNIQAblation.m` |
| `runPNNNPhaseExperiment.m` | `DPD/experiments/runPNNNPhaseExperiment.m` |

Los algoritmos y diagnósticos auxiliares permanecen compartidos en `toolbox/`.
`runPA()` inicia un barrido PA completo y la comparación; `comparePA(P,sweep)`
reutiliza predicciones, reajusta fixed Ridge y exporta. **`comparePA(P)` sin sweep
inicia un barrido completo**, igual que antes. No es una consulta de solo lectura.
Los PA legacy sin `experiment.mat` siguen accesibles con `load`/`readtable`, pero
no admiten reanudación/comparación con el formato actual; no se han migrado.

La copia anterior a esta reorganización es
`C:\Sergi\Investigacion\Códigos\DPD_desarrollo_interno\reorganize_20261014_PA_DPD`.
Incluye el árbol de trabajo completo (también medidas, resultados y WIP), sin el
directorio Git. Su subcarpeta `validation/` conserva las comprobaciones de esta etapa.
La copia de la simplificación anterior y sus tests siguen intactos.

## 1. Qué hace cada flujo

| Flujo | Qué aprende o hace | Entrada principal |
|---|---|---|
| Modelado del PA | Aproxima entrada del PA → salida medida; no genera por sí solo un DPD validado | `runPASweep`, configuración por defecto `experiment20260429T134032_forward_xy.mat` |
| DPD OFFLINE por imitación del ILC | Aprende `u → x_ILC`: referencia deseada → predistorsión calculada por ILC | `trainDPD(measurementFile,P)` |
| Barrido DPD | Compara los tres modelos para varios presupuestos P sobre una captura | `runDPDSweep(measurementFile,20:10:500)` |
| ONLINE | Carga modelos explícitos y genera tres señales, sin ajustar nada | `runDPDOnline(modelFile,inputFile)` → `applyDPDModels` |
| Experimentos auxiliares | Diagnostican o comparan variantes; no sustituyen automáticamente modelos principales | Ejecutores PNNN/PN-IQ separados, §5 |

En la captura DPD real `measurements/experiment20260429T134032_xy.mat`, **x es u e y es x_ILC**, no la salida RF del PA. ONLINE solo necesita x. “ONLINE” significa inferencia de un registro completo, no procesamiento en tiempo real ni transmisión RF.

El preprocesado elimina la media compleja del registro completo cuando `removeDC=true`. La memoria es **periódica**: al principio del registro se usan muestras de su final. No se recortan los transitorios como en un flujo causal continuo. ONLINE aplica el mismo tratamiento a x; no restaura DC de salida, no hace CFR, clipping ni ajuste automático de potencia.

## 2. Datos, ajuste y complejidad

`config/getComparisonConfig.m` contiene la configuración científica. `buildCommonComparisonSplit` selecciona determinísticamente aproximadamente el 10% para identificación; dentro de este conjunto separa 85%/15% para selección interna PNNN, equilibrando amplitudes. La selección histórica usa información del objetivo: no es un muestreo aleatorio independiente de y.

En el barrido DPD real hay 491520 muestras: 49151 de identificación y 442369 en el complemento. Los modelos definitivos se ajustan sobre identificación. El complemento sirve para desarrollo/selección de P; **no es un test externo**, y la NMSE de señal completa incluye identificación.

| Modelo | Entrenamiento | Qué significa P |
|---|---|---|
| GMP-DOMP | Diccionario de 673 regresores complejos; DOMP selecciona soporte y LS ajusta coeficientes | P/2 coeficientes complejos = P escalares reales |
| PN-IQ-GMP | Rota entrada/objetivo a fase de referencia, factoriza GMP, separa características I/Q y elimina ceros estructurales; DOMP selecciona y LS ajusta dos salidas reales | P/2 características compartidas, cada una con coeficiente I y Q |
| PNNN | Características con fase normalizada, red densa, poda y fine-tuning | Pesos activos de ambas capas + 14 sesgos protegidos |

Los lineales principales usan LS sin Ridge sobre columnas escaladas por su pico y objetivo de pico unitario; después se deshace el escalado de coeficientes para predecir. PN-IQ dispone de 1285 candidatos reales tras eliminar 61 componentes Q estructuralmente nulas. En un barrido, cada familia calcula una ruta DOMP y ajusta sus prefijos. Las tres lambdas fixed Ridge son **referencias suplementarias** sobre esos soportes, no los tres modelos del bundle ONLINE.

Los FLOPs son un calendario analítico de inferencia, no tiempos de MATLAB: suma/producto real = 1, producto complejo = 6, suma compleja = 2. Activaciones, exponenciales, raíces y divisiones se informan aparte. Para PNNN se publica el coste disperso ideal, que exige omitir pesos cero en la implementación; `dlnetwork` no garantiza esa aceleración. P310 neuronal cuenta 296 pesos + 14 sesgos y 780 FLOPs de núcleo por muestra según esta convención.

`MaxAbsRealParameter` es un nombre histórico: actualmente GMP usa el **módulo complejo** del coeficiente de comparación, PN-IQ el máximo absoluto I/Q, ambos en coordenadas de pico; PNNN considera pesos/sesgos activos de fc1 y fcOut en coordenadas z-score. No son ganancias físicas directamente comparables. No se ha cambiado su definición en esta auditoría.

## 3. PNNN, paso a paso

1. **Fase:** `r(n)=conj(x(n))/abs(x(n))`, o 1 si x=0. Se rotan los 14 taps `x(n)…x(n-13)` y el objetivo con r(n); al final se restaura la fase multiplicando por `conj(r(n))`.
2. **84 características:** 14 partes reales + 14 imaginarias de taps rotados + 14 envolventes para cada potencia `[1 3 5 7]` (56). El orden de canales forma parte del contrato de los modelos guardados.
3. **Red:** `84 → fc1(12) → sigmoide → fcOut(2)`. Hay 1032 pesos y 14 sesgos: 1046 parámetros densos. Las dos salidas representan I/Q del objetivo rotado. Entradas y salidas se normalizan por z-score con estadísticas guardadas.
4. **Elegir épocas:** `trainDensePNNN` entrena una red provisional con entrenamiento/validación internos. Escoge la época de menor pérdida de validación; selecciona también una duración de fine-tuning usando el presupuesto de referencia P200. Los máximos se escalan para conservar el presupuesto histórico de actualizaciones: 150 en configuración no significa necesariamente 150 épocas del subconjunto actual.
5. **Fuente definitiva:** se descarta la red provisional, se reinicializa con la misma semilla y se recalculan estadísticas sobre toda identificación. Se entrena durante las épocas elegidas. Esta fuente densa se guarda y no se modifica al recorrer P.
6. **Poda:** `pruneAndFineTunePNNN` ordena globalmente por magnitud los pesos de ambas capas y deja P−14 activos. Protege todos los sesgos. Las máscaras mantienen a cero pesos y gradientes podados; se reaplican después de Adam para impedir su reaparición. No se elimina una neurona completa por defecto.
7. **Fine-tuning:** cada P empieza desde la **misma fuente densa**, no desde el P anterior, y usa la duración compartida seleccionada en P200. `fitSparsePNNN` guarda red, máscaras, estadísticas, características y métricas. Reanudar la misma ejecución reutiliza checkpoints; una captura nueva no hereda automáticamente estos elementos.

### Canal 15: histórica frente a corregida

El canal 15 es `Im(r(n)*x(n))`: matemáticamente cero, pero la aritmética compleja deja residuos muy pequeños. El normalizador histórico sustituye sigma por 1 **solo si es exactamente cero**; dividir por una sigma casi nula convierte redondeo en una entrada de orden unidad.

Se investigaron características, pesos/máscaras, rotaciones e intervención con pesos congelados. En P310 histórico las conexiones de ese canal ya estaban podadas: anularlo no cambiaba la predicción. Sí afectaba directamente a la red densa y, menos, a P500. Esto no explicaba por sí solo la brecha de varios dB frente a PN-IQ.

La variante experimental usa `zeroCurrentImaginaryTap=true`: fija ese canal a cero antes del z-score, conserva las otras 83 posiciones y entrena su propia red/estadísticas. La opción ausente o falsa preserva la histórica. No se debe activar sobre un modelo antiguo y presentarlo como modelo reentrenado equivalente.

**El experimento real corregido P310 sí está ejecutado**, en `results/pnnn_phase_experiment/tp1bfa4378_d472_497b_87d7_b5c638c34392/`. Su `comparison.csv` da:

| PNNN P310 | NMSE evaluación | Épocas densa / fine-tuning | Parámetros / FLOPs |
|---|---:|---:|---:|
| Histórica | -35.482326 dB | 1228 / 164 | 310 / 780 |
| Cero estructural, entrenamiento nuevo | -35.525255 dB | 1221 / 163 | 310 / 780 |

La mejora es **0.042929 dB**, no una solución demostrada a la brecha ni una mejora RF. Se mantuvo el criterio de elegir épocas, no su número fijo. Además, un peso contado como activo sobre una entrada cero puede ser funcionalmente inefectivo: P es el presupuesto de máscaras, no necesariamente P grados de libertad efectivos. No se alteraron la poda ni los FLOPs para ocultarlo.

## 4. Qué se reutiliza y qué significa “igual”

Una refactorización debe conservar señales muestra a muestra, coeficientes, normalizaciones, máscaras, soportes, NMSE y conteos. Guardar/recargar o cambiar la organización del código **no** es una nueva comparación científica. Los tests de recarga y las comparaciones históricas usan igualdad exacta cuando corresponde; las tolerancias de ciertas pruebas sintéticas de rotación reflejan la precisión simple de la red.

En cambio, IQ sin PN, soporte diferente o una red entrenada con cero estructural son **experimentos distintos**: se espera que puedan cambiar las predicciones. El control consiste en mantener captura, split, presupuesto y protocolo pertinentes, no en exigir igualdad de NMSE. Una mejor equivariancia tampoco implica por sí sola menor NMSE o mejor ACLR.

El selector requiere todos los puntos cada 10 parámetros dentro de la ventana hacia delante de 100, incluido su extremo, y mejora máxima conjunta ≤0.20 dB. El piloto `[20 100 200 340]` daba un falso P200 por falta de cobertura; ahora se conserva el piloto pero no se selecciona sin evidencia. El histórico PA sigue seleccionando **P340**; el barrido DPD completo selecciona **P310**. Son problemas y métricas distintos.

## 5. Ablación PN-IQ y mapa de código

Se compararon cinco variantes en P100/P310/P500: GMP histórico, PN-IQ histórico, IQ nativo, IQ restringido a candidatos comunes y PN-IQ ajustado sobre el soporte elegido por ese IQ común. “Mismo soporte” significa los mismos pares `(SourceRegressorIndex, IsQ)`, **no los mismos valores de las columnas** tras cambiar la representación. No se trasladaron forzosamente características estructuralmente nulas.

PN-IQ ganó **0.4952 / 0.4481 / 0.5970 dB** sobre IQ con soporte común. Frente a GMP histórico perdió 0.3381 dB en P100 y ganó 0.4606/0.5358 dB en P310/P500. Hay interacciones de representación, libertad I/Q y selección: no es una mejora atribuible exclusivamente a “rotar la fase”. La SVD mostró rango completo y condiciones hasta `3.59e8`; la perturbación dirigida del objetivo cambió NMSE como máximo `1.83e-10 dB`. Es evidencia numérica limitada en doble precisión, no garantía universal.

| Clase | Archivos/carpetas | Decisión de auditoría |
|---|---|---|
| Principal DPD | `trainDPD`, `runDPDSweep`, `runDPDOnline`, `applyDPDModels`, `exampleDPD.m` | Mantener interfaces y bundles |
| Ciencia compartida | `toolbox/domp`, `splits`, `pn_gmp_comparison`, `pnnn`; ajustes/contadores en `toolbox/sweep` y `complexity` | No cambiar orden de operaciones, semillas ni contratos |
| PA/presentación | `runPASweep`, `comparePA`, exportadores y `plotDPDSweepFigures` | Mantener; no confundir con entrenamiento DPD ni con un mero visor sin efectos |
| Experimental PNNN | `runPNNNPhaseExperiment`, `diagnosePNNNPhaseNormalization`, `auditPNNNPhaseSource` | Conservar: entrenamiento explícito, diagnóstico congelado y auditoría de fuente/máscaras son tareas distintas |
| Experimental PN-IQ | `runPNIQAblation`, `predictPNIQAblation`, `diagnosePNIQAblation` | Conservar: el predictor IQ/auxiliares no es un duplicado sustituible por ONLINE de tres modelos |
| Referencia separada | `pn-iq-gmp/pn-iq-gmp/` | Implementación autocontenida no referenciada por el flujo principal revisado; no borrarla ni añadirla recursivamente al path principal |
| Histórico científico | `results/`, informes y checkpoints | Conservados íntegramente; necesarios para reproducir/comparar |
| Desarrollo interno | Tests y fixtures sintéticos | Copia recuperable externa antes de retirarlos del árbol público; véase §7 |
| Residuo accidental | Archivo raíz `= @(git diff --cached --name-only)` | Ayuda de `less`, sin referencias en el pipeline; retirado con copia externa |

### Limpieza aplicada y lo que no se ha “simplificado”

- `runDPDOnline` carga solo x del archivo de entrada; no materializa y ni variables ajenas que no utiliza.
- `diagnosePNNNPhaseNormalization` reutiliza la primera inferencia, cuyo ángulo obligatorio es 0°, como baseline. Elimina una llamada completa redundante sin cambiar tablas ni señales.
- Ridge calcula los máximos directamente al ajustar, sin tres matrices auxiliares ni un segundo recorrido. Parámetros/FLOPs se consultan directamente en las tablas de los soportes originales.
- Ambos barridos solo preparan características y fuente densa PNNN cuando falta algún punto. PA y DPD conservan coordinadores separados.
- PNNN mantiene comprobaciones de presupuesto, sesgos y pesos podados, pero no revalida en cada minibatch tipos y máscaras construidos internamente. La construcción de características queda en una sola función, con el mismo orden de operaciones.
- La población GMP se canoniza una vez; se retiran su salida de depuración no usada y el script de cuatro constantes `modelsConfig.m`, ahora explícitas en el constructor. Se conservan orden y definiciones de candidatos.
- Se eliminan helpers de una operación, alias, conversiones/recuentos repetidos y comentarios redundantes. En SVD ya no se resuelve dos veces el mismo LS solo para comparar el resultado consigo mismo; se conservan refit, perturbación y comprobación contra la predicción guardada.

Se conservaron deliberadamente: la compactación/remapeo de columnas Q en bundles de diferentes P; las máscaras antes/después de Adam; validaciones de presupuesto y protección de archivos; el tratamiento por familia para no convertir predicciones lineales double a single al concatenar PNNN; y el guardado antes de figuras/diagnósticos, que permite recuperar resultados si estos fallan. No son duplicaciones prescindibles.

PA y DPD comparten ajustes, pero sus coordinadores conservan protocolos de evaluación y exportación distintos. No se vectorizaron los bucles científicos ni se cambió el normalizador para reducir líneas. La referencia `pn-iq-gmp/pn-iq-gmp/` permanece separada; se propone excluirla de la distribución pública, no borrarla del archivo científico local.

### Documentación histórica que no debe tomarse como estado actual

El plan, la guía ILC y el informe PNNN tienen notas de actualización: los 49 bundles DPD y el experimento corregido ya existen. Las secciones originales conservan la cronología, no son tareas que deban repetirse. `testDPD` aparece en narraciones históricas, pero **ya no existe como archivo en el árbol actual**: no es un comando operativo ni se ha recreado.

## 6. Comandos que importan

Ejecutar desde la raíz indicada tras `addpath('PA','DPD','PA/experiments','DPD/experiments')`. No usar `addpath(genpath(pwd))`: incluiría la referencia autocontenida y podría ocultar funciones homónimas.

Para una nueva captura con la misma convención física, tras decidir P y autorizar su entrenamiento:

```matlab
modelFile = trainDPD('measurements/experiment20260429T134032_xy.mat',310);
[outputFile,dpd] = runDPDOnline(modelFile,'measurements/experiment20260429T134032_xy.mat');
```

Sustituir las rutas por los nuevos MAT. El primero inicia un entrenamiento completo en una carpeta nueva; el segundo solo genera señales. `dpd(1:3).yvalmod` contiene GMP, PN-IQ y PNNN, con sus etiquetas `modeltype`. `exampleDPD.m` inicia OFFLINE y ONLINE con P310; para reutilizar el modelo guardado, comentar `trainDPD` y activar la ruta indicada en el script.

Para usar ahora el P310 histórico, sin entrenar:

```matlab
runDir = 'results/dpd_sweep/tp98e6f7ac_f3cd_44c7_9ce3_fd8619319748';
modelFile = fullfile(runDir,'models_P0310.mat');
[outputFile,dpd] = runDPDOnline(modelFile,'measurements/experiment20260429T134032_xy.mat');
```

Para cargar/inspeccionar e inferir en memoria:

```matlab
addpath(genpath('toolbox'));
S = load(modelFile,'bundle');
X = load('measurements/experiment20260429T134032_xy.mat','x');
dpd = applyDPDModels(S.bundle,X.x);  % No escribe archivos sin tercer argumento.
```

Barrido nuevo y reanudación son llamadas diferentes; no ejecutarlas solo para consultar:

```matlab
% NUEVO: costoso; entrena una fuente densa y ajusta todos los presupuestos.
sweep = runDPDSweep('measurements/experiment20260429T134032_xy.mat',20:10:500);
% REANUDAR la misma ejecución: usa sus checkpoints y completa lo que falte.
sweep = runDPDSweep(runDir);
% CONSULTAR: no ajusta, exporta ni modifica nada.
S = load(fullfile(runDir,'dpd_complexity_sweep.mat'),'sweep');
S.sweep.results
S.sweep.selection
```

`experiment.mat` fija datos/configuración/split; `linear_sweep.mat` guarda rutas y modelos lineales; `sweep_dense_source.mat` la densa común; `pnnn_target_XXXX.mat` cada poda; `models_PXXXX.mat` el bundle autocontenido ONLINE. `rf_measurements.csv` relaciona posteriormente cada modelo con su medida RF y se preserva al reanudar. La reanudación puede regenerar tablas/figuras: para inspeccionar, usar `load`.

El P340 recuperado sigue en `results/dpd_offline/P340_20261011_210433_267_ed6696d7-1826-412c-b689-f7f61a53abfa/models.mat`. El bundle corregido está en la carpeta experimental del §3; cargarlo debe ser una decisión explícita, no una sustitución silenciosa.

## 7. Validación y recuperación del desarrollo interno

Se conservó una instantánea anterior a esta refactorización en:

`C:\Sergi\Investigacion\Códigos\DPD_desarrollo_interno\refactor_20261013_221701`

Incluye código, documentación y todos los tests originales, también los no versionados y el WIP. `tests_validated/` conserva además las pruebas ajustadas al refactor: los 24 archivos se verificaron byte a byte antes de retirar `tests/`. Para volver a ejecutarlas, copiar esa carpeta como `tests/` en la raíz del repositorio; esas pruebas corresponden al estado anterior a la reorganización PA/DPD y requieren actualizar nombres y rutas según el mapeo siguiente. Los ejecutores científicos y sus informes no se retiran.

Las comparaciones usan igualdad exacta, no una tolerancia nueva para aceptar cambios del refactor:

| Comprobación | Resultado |
|---|---|
| Regresión sintética antes/después | Iguales población GMP, modelos, coeficientes, Ridge, predicciones, tablas, red densa, normalización, épocas y máscaras |
| DPD histórico P310/P340, 491520 muestras | Seis señales idénticas; NMSE, métricas guardadas, presupuestos, metadatos, soportes y máscaras iguales |
| Ablación P100/P310/P500 | Los 15 predictores reproducen exactamente las señales guardadas; sin repetir DOMP |
| SVD y sensibilidad, fixture P12 | Tablas antes/después idénticas; no se repitió la SVD de las capturas reales |
| Selección | DPD P310, PA P340; piloto disperso rechazado |
| PA P340 guardado | NMSE completa, presupuestos, máscaras y FLOPs comprobados sin ajustar modelos |
| DPD sintético P20/30/40 | Capturas independientes, poda desde fuente común, reanudación completa/parcial y recarga aprobadas |
| PNNN cero estructural y ONLINE en otro MATLAB | Prueba sintética aprobada; recarga histórica y experimental con predicciones idénticas, sin entrenamiento al inferir |
| Ajustes, Ridge, espectros y figuras | Pruebas existentes aprobadas, incluido FIG visible y contratos de exportación |
| `checkcode` | 53 archivos MATLAB propios sin avisos |
| `git diff --check` | Sin errores |

Matiz PA histórico: la NMSE de identificación GMP guardada es −38.903441250607273 dB; recalcularla usando el subconjunto de la predicción completa guardada da −38.903441250607464 dB (1.92e−13 dB). Son evaluaciones históricas distintas —matriz de identificación frente a predicción por bloques—, no una diferencia introducida al refactorizar: se confirmó igualdad exacta usando también los helpers anteriores. No se cambió ninguno de esos valores ni se reajustó el modelo. El checkpoint neuronal PA antiguo guarda predicción y máscara, no la red final: su comprobación es de referencias guardadas, no una nueva inferencia del modelo ausente.

Los logs `dpd_final_*.log`/`.err` y `final_checks.log`/`.err` se conservan junto a la copia externa. Un primer test detectó una lista de rutas cell en la nueva reanudación; se corrigió a string y se repitió la prueba. No se hizo entrenamiento real, barrido científico nuevo ni transmisión RF.

### Balance de la simplificación

Comparación con la instantánea inicial del WIP, no con HEAD. Son líneas físicas MATLAB, incluidos comentarios y blancos; se excluyen `third_party` y la referencia independiente `pn-iq-gmp`.

| Categoría | Archivos antes → después | Líneas antes → después |
|---|---:|---:|
| Código principal y presentación | 47 → 46 | 4476 → 4172 |
| Experimentos y diagnósticos científicos | 7 → 7 | 636 → 635 |
| Tests internos, conservados externamente | 24 → 0 | 2258 → 0 |
| Total MATLAB propio | 78 → 53 | 7370 → 4807 |

La reducción del código conservado es **305 líneas (5.97%)**; otras 2258 líneas se retiraron del árbol principal, pero siguen recuperables. No se presenta esa retirada como una simplificación del algoritmo. Fuera de los tests solo se retiraron `modelsConfig.m` y el archivo accidental de ayuda de `less`. README y documentación se contabilizan aparte y no forman parte de esta reducción MATLAB.

## 8. Qué falta en el laboratorio

Confirmar con el tutor, por captura, semántica de x/y, frecuencia de muestreo, amplitudes/unidades, DC, alineación y periodicidad. Verificar longitud, PAPR/CFR, límites DAC y sincronización antes de transmitir. Elegir explícitamente el bundle histórico o experimental y registrar su procedencia junto con la captura RF.

Después, medir la salida del PA con condiciones comparables para los tres DPD y el baseline pertinente: NMSE RF, espectro/ACLR-ACPR, EVM cuando corresponda y régimen de potencia/back-off. Los FLOPs ideales no sustituyen medidas de latencia/recursos de implementación. La confirmación científica independiente requiere otra captura/campaña con selección y evaluación separadas. **En esta auditoría no se transmitió RF, no se entrenó ninguna red real ni se ejecutó un barrido completo.**

Para detalles, conservar y consultar los informes [PNNN](PNNN_PHASE_NORMALIZATION_DIAGNOSIS_2026-10-12.md), [ablación PN-IQ](PNIQ_GMP_ABLATION_2026-10-12.md) y [flujo DPD](DPD_ILC.md). No hubo commit ni push.

## 9. Auditoría de nombres MATLAB (14/10/2026)

Se revisaron los **53 archivos propios** completos: **49 funciones principales,
51 funciones locales y 7 métodos** (107 declaraciones), además de los dos
scripts y las dos clases. Se excluyen terceros y la referencia independiente
`pn-iq-gmp/`. No se modificaron algoritmos ni campos de datos/modelos.

| Nombre anterior | Nombre actual | Motivo |
|---|---|---|
| `run_linear_sweep` | `fitLinearModels` | Ajusta las dos familias para uno o varios presupuestos; no ejecuta un experimento. |
| `run_fixed_ridge_sweep` | `fitFixedRidge` | Reajusta soportes dados con lambdas fijas. |
| `fit_complex_gmp_domp` | `fitComplexGMP` | Ajusta GMP con selección DOMP o soporte suministrado; el ajuste no siempre ejecuta DOMP. |
| `fit_pniq_gmp` | `fitPNIQGMP` | Identifica la familia y mantiene la convención del auxiliar de ajuste GMP. |
| `fit_sparse_pnnn_target` | `fitSparsePNNN` | Poda y ajusta un presupuesto desde la fuente densa. |
| `prepare_pnnn_dense_source` | `trainDensePNNN` | Selecciona épocas y entrena la red definitiva; no es mera preparación de datos. |
| `GMP_createRegressorManager` | `buildGMPManager` | Construye la población y su gestor, no evalúa una matriz. |
| `sel_indices` | `selectPeakSegment` | Selecciona el segmento histórico usando el pico del objetivo. |
| `selectXYByMapping` | `mapModelSignals` | Asigna las señales según xy/yx; no selecciona muestras. |
| `predictPhaseNorm` | `predictPNNN` | Predice la salida compleja y restaura escala/fase; no devuelve una fase normalizada. |
| `getFairDOMPComparisonConfig` | `getComparisonConfig` | Configura la comparación compartida, incluida PNNN, no únicamente DOMP. |
| `regressorManager.initialization` | `regressorManager.buildPopulation` | Construye candidatos GMP. |
| `regressorManager.removerepeated` | `regressorManager.finalizePopulation` | Añade auxiliares, canoniza, filtra límites y elimina duplicados. |
| `Regressor.sortindexes` | `Regressor.sortLags` | Ordena retardos, no índices de muestras. |
| `fase_pmpi` (local) | `wrapPhase` | Envuelve fase en radianes a [-pi,pi). |
| `iqOperations` (local) | `countIQOperations` | Cuenta operaciones del control IQ. |
| `ridgeFit` (local) | `fitRidge` | Verbo de ajuste coherente con los demás ajustes. |

También se aclaran variables de alcance local: resultado `sweep → linear`
en `fitLinearModels`; `rManagerGMP → manager`; objetivo, fracción e índices
de `selectPeakSegment`; `network`, `normalizedFeatures`,
`normalization` y `phaseRotation` en `predictPNNN`; e índices
`featureIndices` para no confundirlos con la matriz de características.
Los campos guardados y nombres de checkpoints conservan sus nombres originales.

Se mantienen los ejecutores PA/DPD, los experimentos y sus nombres.
También `countModelOperations` (calendarios de operaciones),
`countModelFLOPs` (conversión a FLOPs), `countPNNNFLOPs`
(costes de arquitectura densa) y `countSparsePNNNFLOPs` (coste disperso
ideal): sus contratos ya distinguen las responsabilidades.
`countPNNNParameters` cuenta la arquitectura densa;
`summarizeTrainableParameters` inspecciona pesos/sesgos y máscaras reales.
`buildCommonComparisonSplit` y `compactGMPRegressors` describen
correctamente la partición común y la serialización compacta, respectivamente.
Los auxiliares de figuras conservan sus nombres específicos: generan/exportan
las mismas figuras, sin cambiar datos ni presentación.

**Compatibilidad:** se mantienen `Regressor` y `regressorManager`, sus
constructores y propiedades para no cambiar identidades de clases históricas.
Los bundles actuales guardan regresores compactos; no contienen callbacks
hacia los métodos renombrados. La única llamada dinámica del código propio
es `dlfeval(@modelLoss,...)`, que no cambia. No hay alias ni wrappers antiguos.

La instantánea previa recuperable está en
`C:\Sergi\Investigacion\Códigos\DPD_desarrollo_interno\names_20261014`.
Incluye todo el WIP, resultados y medidas, sin Git. Sus 53 fuentes MATLAB
se verificaron por hash antes de editar. `validation/inventory.md` enumera
todos los archivos y las 107 declaraciones revisadas; `validation/`
conserva las pruebas externas adaptadas y los logs. Las copias anteriores
permanecen intactas.

Las llamadas compartidas quedan así (añadir `config/` y `toolbox/` al
path si se usan sin un ejecutor):
```matlab
cfg = getComparisonConfig(projectRoot);
[linear, models] = fitLinearModels(x,y,split,cfg);
fixed = fitFixedRidge(x,y,split,cfg,linear);
[denseSource,features,targets,rotation] = trainDensePNNN(x,y,split,cfg,P);
point = fitSparsePNNN(denseSource,P,features,targets,rotation,y,split,cfg);
```
`trainDensePNNN` inicia entrenamiento: no debe invocarse solo para consultar
un checkpoint. Las llamadas principales de PA y DPD del README no cambian.

### Verificación de esta auditoría nominal

- MATLAB R2026a, `validation2.log`: **ALL PASS**, salida 0. Los tests permanecen externos.
- Igualdad exacta de las seis señales DPD P310/P340, metadatos, NMSE,
  modelos lineales, pesos, máscaras y presupuestos neuronales.
- Igualdad exacta de las quince predicciones de ablación P100/P310/P500.
- Selección histórica PA P340 y DPD P310 conservada.
- Fixture de 160 muestras, dos épocas: población, ajustes lineales, Ridge,
  predicciones, tablas, normalización, pesos densos, épocas y máscaras idénticos.
- Reanudaciones sintéticas PA/DPD con métricas idénticas y checkpoints existentes;
  `comparePA` desde fuera del repositorio con path limpio; `runPA` automático/manual aprobado.
- Los 53 MATLAB pasan `checkcode`, declaración/archivo y `which` sin colisiones.
  No quedan llamadas operativas a los nombres sustituidos. `git diff --check` pasa.
- 548 archivos ajenos a código/documentación comprobados por SHA-256 sin cambios,
  incluidos medidas, resultados y checkpoints. Los cuerpos de los 53 MATLAB
  coinciden con el estado previo tras aplicar únicamente las sustituciones nominales declaradas.

Balance exclusivo de esta etapa: 26 MATLAB tocados (11 renombrados), tres documentos
actualizados, **53 archivos / 4802 líneas MATLAB antes y después**. No se atribuye
a esta auditoría el WIP de reorganizaciones anteriores. La documentación se cuenta
aparte. `exampleDPD` no se modificó: se corrigió la guía porque su WIP actual inicia
OFFLINE, no una simple recarga.

La primera comprobación externa falló solo por usar una expresión de recuento
incompatible con `regexp` de MATLAB; se corrigió en el comprobador, antes de ejecutar
la regresión numérica. En PA aparecieron los avisos esperados de parametrización
y de ausencia de `info_signal` en el fixture que sí contiene `fs`; no impidieron
la exportación. Sin entrenamiento real, barrido completo, transmisión RF, commit ni push.

