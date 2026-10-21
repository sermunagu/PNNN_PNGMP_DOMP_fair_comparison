# Refactor de simplicidad DPD — 12/10/2026

> Organización actual (14/10/2026): ejecutores PA en `PA/` y DPD en `DPD/`;
> `runPNIQAblation` y `runPNNNPhaseExperiment` están en `DPD/experiments/`.
> Antes de usar ejemplos: `addpath('PA','DPD','PA/experiments','DPD/experiments')`.
> Las narraciones históricas conservan sus nombres; el [mapeo actual](GUIA_REPOSITORIO_DPD.md#entradas-reorganizadas-14102026)
> indica sus sustitutos, sin cambiar medidas, modelos ni resultados.

## A. Alcance y cambios

Repositorio original: `C:\Sergi\Investigacion\Códigos\PNNN_PNGMP_DOMP_fair_comparison`.
Rama `main`, HEAD `20eaf4cee23b713a47d24fe92bc4159ab63b597e`.
Se leyó íntegramente la auditoría de simplicidad y se inspeccionaron fuentes y
consumidores. No se hicieron commits ni pushes.

WIP previo conservado: eliminación de `testDPD.m`, auditoría no versionada y
carpeta `pn-iq-gmp/`. No se modificaron medidas, MAT históricos, modelos P340,
resultados versionados ni el submódulo matlab2tikz.

- Barrido: carpeta nueva por ejecución; reanudación explícita por carpeta.
  `experiment.mat` conserva cfg, x/y preprocesadas y split; cada etapa/presupuesto
  tiene su MAT directo. No hay identidad criptográfica, schema ni migración automática.
- Se guardan resultados numéricos antes de exportar. Eliminados preflight gráfico,
  escritura atómica CSV/TXT y copias redundantes de datos de selección.
- Comparación seleccionada: carga directa, predicciones guardadas y reajuste
  Fixed Ridge sobre soportes guardados. Conservadas todas sus figuras.
- PNNN: sin hashes antes/después de podar. Cada presupuesto sigue recibiendo
  la misma red densa; esta propiedad se verifica comparando valores en tests.
- Fixed Ridge reutiliza `buildPNIQFeatures` en identificación e inferencia.
  Una vista local de la población conserva la correspondencia de descriptores.
- Clases: retirados genética, evaluación heredada, cachés de matrices sin consumidores
  y bases FV/CVS/MP/DDR no utilizadas. Conservados los bucles GMP, canonización,
  filtrado y deduplicación originales.
- Eliminados validadores redundantes de formatos y presets; mantenidas condiciones
  matemáticas relevantes, protección de biases y máscaras de gradientes/pesos.
- La regresión de recuentos pasa a un test con soporte histórico explícito,
  en vez de condicionar producción por S/K/F. Se mantiene la aritmética del contador.
- Probe WL-GMP: ahora recibe `referenceFile` explícito; no se ejecutó el experimento.
  El script de conteo denso ya no imprime un supuesto P340 fijo.
- OFFLINE/ONLINE y `llamadas_DPD.m` no se modificaron.

Archivos eliminados: `buildExperimentSignature.m`, `buildNetworkSignature.m`,
`preflightPaperFigureToolchain.m` y el test exclusivo de firmas.

Funciones retiradas, agrupadas:

- Barrido: `buildSweepIdentity`, `resolveResultDirectory`,
  `writeTableAtomically`, `writeTextAtomically`, `deleteIfPresent`,
  `hasCurrentLinearProtocol`, `hasCurrentFixedProtocol`.
- Firmas: las dos funciones públicas y sus auxiliares de serialización,
  endianess, actualización/finalización del digest; preflight y `removeDirectory`.
- Regressor: `mutatememory`, `mutateordershuffle`, `mutateorderincrement`,
  `crossover`, `crossoverkk`, `buildRegressor`.
- Manager: `buildX`, `printModel`, `evaluation`, `selection`, `crossover`,
  `mutation`, `buildUcustomX`, `clearRegressors`, `regress`, `prepareForSave`.
- Validadores/helpers: `validateConvention`, `validateSupport`,
  `validateInputSignal`, `validateOrders`, `periodicTapIndices` (inline),
  `applyDefaults`. La selección de columnas canónicas sigue existiendo como
  `canonicalSpecifications`, sin su validación de schema.

Dependencias retiradas del flujo: Java/SHA-256 y las llamadas sin consumidores
a evaluadores/bases genéticas externas. matlab2tikz sigue siendo dependencia
de exportación, no de identificación.

## B. Líneas antes/después

Recuento de líneas físicas incluyendo comentarios y blancos internos, sin
contar el salto vacío final. Solo cambios de esta tarea; excluye el WIP previo.

| Archivo | Antes | Después |
|---|---:|---:|
| `config/getFairDOMPComparisonConfig.m` | 93 | 91 |
| `count_dense_parameters.m` | 84 | 81 |
| `docs/DPD_ILC.md` | 183 | 139 |
| `docs/PLAN_DPD_TRANSACTIONS_2027.md` | 344 | 324 |
| `run_parameter_sweep.m` | 505 | 167 |
| `run_selected_comparison.m` | 524 | 396 |
| `run_widely_linear_gmp_probe.m` | 99 | 98 |
| `tests/run_dpd_model_reuse_test.m` | 101 | 101 |
| `tests/run_experiment_signature_test.m` | 47 | 0 |
| `tests/run_export_paper_figure_roundtrip_test.m` | 67 | 64 |
| `tests/run_operating_point_selection_test.m` | 94 | 94 |
| `tests/run_pnnn_shared_dense_sweep_test.m` | 170 | 169 |
| `tests/run_scientific_pipeline_smoke_test.m` | 65 | 78 |
| `tests/run_simplification_regression_test.m` | 0 | 52 |
| `tests/run_sweep_configuration_test.m` | 53 | 52 |
| `tests/run_sweep_resume_test.m` | 0 | 54 |
| `toolbox/complexity/countModelFLOPs.m` | 179 | 115 |
| `toolbox/pn_gmp_comparison/GMP_createRegressorManager.m` | 27 | 9 |
| `toolbox/pn_gmp_comparison/Regressor.m` | 286 | 100 |
| `toolbox/pn_gmp_comparison/countModelOperations.m` | 258 | 224 |
| `toolbox/pn_gmp_comparison/factorizeGMPRegressor.m` | 92 | 87 |
| `toolbox/pn_gmp_comparison/modelsConfig.m` | 30 | 4 |
| `toolbox/pn_gmp_comparison/regressorManager.m` | 436 | 135 |
| `toolbox/pnnn/buildExperimentSignature.m` | 97 | 0 |
| `toolbox/pnnn/buildNetworkSignature.m` | 62 | 0 |
| `toolbox/pnnn/buildPhaseNormDataset.m` | 128 | 89 |
| `toolbox/pnnn/fit_sparse_pnnn_target.m` | 87 | 74 |
| `toolbox/pnnn/prepare_pnnn_dense_source.m` | 244 | 234 |
| `toolbox/sweep/computeSelectedPointSpectra.m` | 58 | 40 |
| `toolbox/sweep/exportPaperFigure.m` | 60 | 47 |
| `toolbox/sweep/plotSweepPaperFigure.m` | 171 | 140 |
| `toolbox/sweep/preflightPaperFigureToolchain.m` | 55 | 0 |
| `toolbox/sweep/run_fixed_ridge_sweep.m` | 354 | 293 |
| `toolbox/sweep/selectOperatingPoint.m` | 217 | 170 |
| `docs/informe_refactor_simplicidad_2026-10-12.md` (este informe) | 0 | 183 |

Producción modificada: **4146 → 2594**, reducción neta **1552 líneas**.
Tests afectados/nuevos: **597 → 664** (+67 netas); documentación operativa:
**527 → 463** (−64). Antes de este informe: **1909 líneas retiradas, 360 añadidas**,
reducción neta **1549**. Incluyendo este informe: 543 añadidas y 1366
de reducción neta. No se atribuyen a esta tarea las 165 líneas de testDPD ya borradas.

Adiciones relevantes: flujo directo de checkpoints (82 líneas nuevas en el diff,
frente a 420 retiradas), carga directa de comparación (24 frente a 152),
dos tests nuevos (52 + 54 líneas), fixture histórico de operaciones (13),
comparaciones de learnables en tests y documentación operativa.
No se añadieron clases ni helpers productivos nuevos para sustituir los retirados.

## C. Verificación ejecutada

MATLAB R2026a, procesos acotados a 60–300 segundos según prueba.
Referencia anterior a cualquier cambio productivo:
`%TEMP%/dpd_before_refactor_20261012.mat`.

**Comparación antes/después: igualdad exacta (`isequaln`)** de población ordenada
de 673 regresores, tablas lineales, rutas DOMP, correspondencias I/Q, coeficientes
lineales, predicciones, tablas/predicciones Fixed Ridge, NMSE, lambdas, parámetros,
FLOPs, pesos/normalización de la densa, épocas seleccionadas, máscaras, pesos
sparse y predicciones PNNN. Fixture de 160 muestras, presupuestos 20 y 30,
PNNN reducida para la prueba. No es una afirmación sobre un barrido real no ejecutado.

**P340 real:** carga del `models.mat` existente y comparación individual con
`online/inference_20261011_212149_467_46617fa4-db5e-477a-8571-1ca7278bcc17/dpd_signals.mat`.
Las tres señales conservan tipo numérico y valores exactos: error relativo **0**.
Perfil de inferencia sin DOMP, trainnet ni poda. No se guardó sobre ningún resultado
real; tamaños y fechas de los MAT originales se conservaron.

Pruebas completadas (19):

- `run_simplification_regression_test`: referencia antes/después exacta.
- `run_sweep_resume_test`: ejecución nueva, reanudación completa sin entrenamiento,
  un checkpoint ausente → solo una poda/fine-tuning; filas/máscaras/predicciones idénticas.
- `run_domp_unit_test`, `run_linear_complexity_sweep_test`,
  `run_fixed_lambda_linear_sweep_test`, `run_max_abs_real_parameter_test`.
- `run_pnnn_shared_dense_sweep_test`, `run_dpd_model_reuse_test`.
- `run_dpd_offline_online_test`: dos capturas independientes, P20 exacto en
  tres modelos; cuatro DOMP y cuatro trainnet en OFFLINE, cero en ONLINE.
- `run_dpd_online_fresh_session_test`: otro proceso MATLAB, misma salida exacta.
- `run_scientific_pipeline_smoke_test`, `run_sweep_configuration_test`,
  `run_sweep_presentation_test`, `run_measurement_sample_rate_test`.
- `run_operating_point_selection_test`, `run_main_sweep_and_comparison_test`.
- `run_selected_point_spectrum_test`, `run_export_paper_figure_roundtrip_test`,
  `run_paper_figure_contract_test`.

Figuras: NMSE/parámetros, NMSE/FLOPs, rango, espectros principales/Ridge,
AM/PM, error de fase y tiempo generados en fixtures temporales. FIG abierto
visible; contratos de colores, límites, frecuencia, etiquetas y FIG/PNG/TikZ pasan.
Los fixtures gráficos se limpian al terminar: no se regeneran resultados históricos.

`checkcode`: 28 archivos revisados; cinco avisos heredados en regressorManager
(tres NBRAK2, dos AGROW), ningún error de sintaxis. No se reformateó esa lógica
solo para silenciar avisos. `git diff --check`: correcto.

Problema encontrado y corregido: el test de selección apuntaba a
`legacy/sweep_d113e389ab78`, ausente. Usa ahora `legacy/sweep_4329deabf42b`,
que conserva la referencia P340 y pasa las mismas aserciones.
La antigua condición S=100/K=17/F=172 no identifica por sí sola un soporte único:
el test de operaciones usa un soporte concreto disponible, no aquella condición
genérica en producción.

## D. Conservado y pendiente

- Intactos DOMP, split, ajustes principales, arquitectura, semillas, entrenamiento,
  selección de épocas, poda, máscaras y fine-tuning.
- Las dos rutas del contador FLOP siguen siendo necesarias: tablas analíticas
  lineales y especificación PNNN. Sus fórmulas no se cambiaron.
- Se conservan firmas algebraicas I/Q de descriptores: deduplican características;
  no son hashes de archivos/redes y sí intervienen en el recuento científico.
- Se conserva fs/info_signal.fsovs: el MAT PA real carece de fs, mientras el DPD
  tiene ambos. No se puede eliminar esa alternativa sin romper el experimento.
- Se conservan condiciones de selección sin punto estabilizado y monotonía FLOP,
  por su significado científico, así como limpieza de figuras/restauración Visible.
- Nuevos checkpoints no reanudan automáticamente formatos antiguos. Para
  reproducir el flujo antiguo se usa su revisión histórica; los datos permanecen intactos.
- No se ejecutaron entrenamientos reales largos, barrido 20:10:500, probe WL-GMP
  ni RF. No se afirma equivalencia de nuevas campañas no evaluadas.
- Sigma PNNN casi nula, holdout independiente, compatibilidad física y evaluación
  RF siguen pendientes: resolverlos cambiaría el protocolo y requiere otra tarea.
