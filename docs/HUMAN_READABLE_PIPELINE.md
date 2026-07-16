# Pipeline MATLAB legible: arquitectura y auditoría

Esta rama reescribe el flujo del sweep sin cambiar el protocolo científico.
El split actual, sus semillas, el 10 % de identificación, las lambdas, DOMP,
PN-DOMP, la PNNN N12, la poda, los FLOPs y los espectros permanecen iguales.
La revisión científica del split queda expresamente pendiente.

## Métricas

| Métrica | Antes | Después |
|---|---:|---:|
| Archivos `.m` del repositorio | 90 | 66 |
| Funciones MATLAB propias | 308 | 221 |
| Líneas MATLAB | 11 217 | 7 710 |
| Archivos del pipeline ejecutable | 12 | 9 |
| Funciones dentro del pipeline | 95 | 68 |
| Líneas del pipeline | 2 536 | 1 944 |
| Profundidad máxima aproximada del pipeline | 4 | 3 |

Las cifras incluyen código y tests, pero excluyen resultados generados.

## Mapa de decisiones

| Área anterior | Decisión | Motivo |
|---|---|---|
| `main_sweep_and_comparison.m` | KEEP | Entrada interactiva mínima. |
| `run_parameter_sweep.m` | MERGE | Integra tabla, CSV y figuras del sweep en orden científico. |
| `run_fair_PNNN_vs_PNGMP_DOMP.m` | MERGE/RENAME | Sustituido por `run_selected_comparison.m`. |
| `runLinearComplexitySweep.m` | RENAME | Algoritmo completo conservado como `run_linear_sweep.m`. |
| `runFixedLambdaLinearSweep.m` | RENAME | Ridge completo conservado como `run_fixed_ridge_sweep.m`. |
| `runPNNNComparisonStudy.m` | SIMPLIFY | Eliminados H4 y el modo histórico; queda la fuente N12 usada. |
| `runPNNNSparseSweep.m` | RENAME | Ajuste sparse completo conservado por target. |
| `buildSelectedFixedLambdaPredictions.m` | INLINE | Era un wrapper de un único target. |
| `writeSweepPresentationOutputs.m` | INLINE | Tabla y figuras se leen al final del sweep. |
| `writeSelectedPointSpectra.m` | INLINE | Los cuatro espectros se leen junto al punto seleccionado. |
| `computeSelectedPointSpectra.m` | KEEP | Operación Welch común, matemática y reutilizable. |
| `updateSweepCheckpoint.m` | KEEP | Escritura atómica e identidad de artefactos externos. |
| Comparación histórica de seis modelos | DELETE | No pertenece a las tres familias actuales. |
| Helpers coupled/no-PN históricos | DELETE | Solo alimentaban la ruta histórica retirada. |
| Tests de wrappers eliminados | DELETE | Sustituidos por tests de contratos científicos. |

## Lectura recomendada

1. `main_sweep_and_comparison.m`: ejecuta el sweep, solicita el presupuesto y
   abre la comparación seleccionada.
2. `run_parameter_sweep.m`: muestra de arriba abajo carga, split, modelos
   lineales, Ridge, PNNN, checkpoints, tablas y figuras.
3. `toolbox/sweep/run_linear_sweep.m`: contiene completos Complex GMP-DOMP y
   PN-IQ PN-DOMP: matrices, soporte, lambda, refit, predicción y métricas.
4. `toolbox/sweep/run_fixed_ridge_sweep.m`: reajusta las dos familias sobre
   los mismos soportes para `1e-3`, `1e-4` y `1e-5`.
5. `toolbox/pnnn/prepare_pnnn_dense_source.m`: selecciona épocas y reajusta una
   única fuente densa N12.
6. `toolbox/pnnn/fit_sparse_pnnn_target.m`: poda, reajusta, predice y cuenta el
   coste de un presupuesto sparse.
7. `run_selected_comparison.m`: carga los tres modelos, reconstruye las seis
   Ridge del target y escribe las cuatro figuras espectrales.

## Operaciones algorítmicas conservadas

- `GMP_createRegressorManager` y `buildGMPRegressorRows` para la población GMP.
- `selectDOMPSupport` y `selectSharedIQFeatures` para DOMP y PN-DOMP.
- `buildPhaseNormalizedIQRegressors`,
  `removeStructurallyZeroQFeatures` y `computePhaseNormGMPRotation`.
- `fitComplexGMPGrid` y `fitIndependentIQGMP`.
- `buildPhaseNormDataset`, `fitFairPNNNDenseValidation`,
  `refitFairPNNNDense` y `refitFairPNNNSparse`.
- Contadores de parámetros/FLOPs y funciones de pruning.

## Contrato científico visible

La ruta PN-IQ muestra en `run_linear_sweep.m`:

1. construcción de la población GMP compartida;
2. creación de features I/Q phase-normalized;
3. eliminación única de columnas Q estructuralmente nulas;
4. PN-DOMP sobre internal train e identificación;
5. ajuste independiente de coeficientes I y Q;
6. reconstrucción compleja normalizada;
7. restauración mediante la fase conjugada de la fila;
8. predicción sobre identificación y señal completa.

## Validación

Los tests finales comprueban split y rotación, DOMP, firmas, sweep lineal,
Ridge, PNNN sparse, checkpoints, comparación seleccionada, cuatro figuras y
Welch común. La ejecución real reutiliza los artefactos de
`results/parameter_sweep/sweep_0c53cac19fae` y reproduce el target 340.
