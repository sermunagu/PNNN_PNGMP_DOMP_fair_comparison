# Auditoría exhaustiva de simplicidad — PNNN_PNGMP_DOMP_fair_comparison

**Repositorio:** https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison  
**Rama/commit auditado:** `main` — `20eaf4cee23b713a47d24fe92bc4159ab63b597e` (12-10-2026 UTC).  
**Objetivo:** encontrar sobreingeniería, verificaciones innecesarias, abstracciones, duplicaciones y residuos de código heredado, **sin confundir complejidad algorítmica real con complejidad artificial**. No se han modificado archivos.

## 1. Alcance, inventario y criterio

Inventario completo por el árbol recursivo de GitHub: **275 archivos versionados**, con **68 archivos `.m`** (50 de aplicación/algoritmos y 18 tests), **4 documentos `.md`**, **108 `.mat` (punteros Git LFS)**, 13 PDF, 22 FIG, 22 PNG, 22 TikZ, 10 CSV, 2 TXT, un TEX y 3 archivos de configuración Git. Hay **195 archivos bajo `results/`** y un submódulo externo en `third_party/matlab2tikz`.

Se inspeccionó individualmente el contenido fuente de los 68 archivos MATLAB mediante GitHub, se examinó la documentación versionada y se inventariaron los demás ficheros. Los binarios de resultados, medidas LFS, PDF, FIG y el código del submódulo externo **no se han auditado internamente**: no contienen código MATLAB propio que permita evaluar su estilo a partir del árbol, y no se dispuso aquí de MATLAB para ejecutar experimentos. La carpeta local no versionada `pn-iq-gmp/` que apareció antes en `git status` tampoco forma parte de esta auditoría de GitHub.

**Criterio:** la lectura de medidas del tutor se resuelve mediante `load` y acceso a `x`, `y` y, si es necesario, `fs`. Si los datos fallan, MATLAB debe mostrar su error normal. No son necesarios comprobadores universales de formato, hashes de cada captura, UUID, gestores de directorios ni defensas frente a entradas imaginarias. En cambio, no se consideran prescindibles la selección DOMP, la construcción PN-IQ, la normalización, los protocolos de entrenamiento, la poda, las máscaras que impiden regenerar pesos y el cómputo definido de FLOPs. Tampoco es razonable eliminar checkpoints de un barrido que consume horas solo para reducir líneas.

### Síntesis

**Resultado: el repositorio todavía contiene sobreingeniería apreciable, concentrada en el antiguo pipeline de barridos y en mecanismos de identidad/reproducción.** El reciente OFFLINE/ONLINE es sencillo en su fachada, pero sigue llamando indirectamente a código de firmas SHA-256 basado en Java a través de la PNNN. También existe duplicación literal importante de la construcción PN-IQ, código heredado de algoritmos genéticos aparentemente no utilizado y un preflight de exportación gráfica que bloquea el entrenamiento.

## 2. Hallazgos específicos, ordenados por impacto

### A1 — Muy alta: `run_parameter_sweep.m` concentra demasiada infraestructura (506 líneas)

[Ver archivo, L1-L250](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/run_parameter_sweep.m#L1-L250) · [L251-L506](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/run_parameter_sweep.m#L251-L506)

- **L25-L28:** antes de cargar las medidas o entrenar, crea una figura ficticia para comprobar la exportación `matlab2tikz`. Un fallo del exportador impide obtener resultados numéricos aunque el entrenamiento sea perfectamente posible. Separar el cálculo de las figuras; no hacer el preflight obligatorio.
- **L36 y L371-L416:** genera primero `buildExperimentSignature` y después una segunda identidad con configuración, `jsonencode` y `java.security.MessageDigest`. Se puede identificar una ejecución por carpeta/configuración sin calcular hashes criptográficos sobre todo el experimento.
- **L55-L196:** reanudación mediante múltiples checkpoints y validación de identidades, nombres, versiones de protocolo y objetivos. **Conservar** el checkpoint por PNNN y la posibilidad de reanudar barridos. **Simplificar** el formato a `load`/`save` de resultados del mismo barrido, sin firmas, jerarquías de protocolo y comprobaciones de cada campo.
- **L219-L259:** verificaciones pormenorizadas sobre número/estructura de filas, lambdas, parámetros y 13/8 columnas. Las invariantes científicas relevantes deberían estar probadas por los tests; no es necesario convertir cada ejecución normal en una batería de pruebas.
- **L327-L356:** el resultado final expone numerosas copias de datos y subcampos derivados de `selection`, además de duplicar datos de checkpoints.
- **L428-L460:** se implementa una transacción de escritura mediante `tempname` + `onCleanup` + `movefile` + `deleteIfPresent`, para CSV y TXT. Para este entorno basta normalmente `writetable`/`writelines`.
- **L462-L505:** dos funciones de comprobación de versiones de protocolo con abundante `isfield` y comparaciones; candidatas claras a eliminación si se abandona la compatibilidad retroactiva automática.

No propondría dejarlo necesariamente en 40 líneas: **el barrido científico y los checkpoints sí hacen trabajo real**. Propondría separar mentalmente cálculo, puntos de guardado y figuras, recortando primero todo el sistema de identidad, compatibilidad y escritura atómica.

### A2 — Muy alta: Java/SHA-256 también subsiste dentro del OFFLINE nuevo

- [`toolbox/pnnn/buildExperimentSignature.m`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/pnnn/buildExperimentSignature.m): **98 líneas**, serializa decenas de opciones de configuración y calcula tres hashes mediante Java, además de utilidades propias de bytes/endianess.
- [`toolbox/pnnn/buildNetworkSignature.m`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/pnnn/buildNetworkSignature.m): **63 líneas** de SHA-256 sobre pesos de red y estadísticas de normalización, con funciones auxiliares casi paralelas a las del anterior.
- [`toolbox/pnnn/prepare_pnnn_dense_source.m#L131-L149`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/pnnn/prepare_pnnn_dense_source.m#L131-L149): llama a `buildNetworkSignature` al crear la fuente densa.
- [`toolbox/pnnn/fit_sparse_pnnn_target.m#L6-L33`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/pnnn/fit_sparse_pnnn_target.m#L6-L33): **calcula de nuevo la firma dos veces**, antes y después de podar, solo para comprobar que no ha cambiado la fuente. Esto **se ejecuta en el OFFLINE nuevo**, aunque `trainDPDOffline.m` solo tenga ~59 líneas incluyendo blancos.

**Propuesta:** si se abandona deliberadamente la detección criptográfica de modificaciones de checkpoints, suprimir las firmas y estas dos comprobaciones. Mantener siempre el contrato matemático de partir de la misma red densa para cada P. Eliminar campos `digest` y sus comparaciones exige adaptar coherentemente checkpoints, comparación seleccionada y tests: no basta borrar una función aislada. No tocar checkpoints históricos.

### A3 — Muy alta: `run_selected_comparison.m` (525 líneas) vuelve a verificar el mismo experimento continuamente

[Ver L18-L175](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/run_selected_comparison.m#L18-L175)

- **L43-L107:** cuatro cargas de checkpoints con comprobación casi idéntica de `sweepIdentity`, seguida de comparación exacta entre datos que ya estaban en memoria y los guardados (`isequaln` para resumen, tablas y estadísticas). Más código de control que de operación científica durante este bloque.
- **L112-L120:** vuelve a hashear la captura completa para contrastarla con el barrido.
- **L131-L173:** repite el ajuste Ridge sobre soportes guardados y compara fila a fila resultados/frecuencias/FLOPs con múltiples tolerancias. Parte del reajuste es legítima porque produce predicciones para graficar; **la auditoría comparativa exhaustiva cada vez que se hace un gráfico no lo es**.
- **L177-L246 y L248-L525:** presentación, tablas, espectros, AM/PM, fase y figuras. Muchas líneas aquí **sí producen las visualizaciones solicitadas para el artículo**. No considerar todas ellas prescindibles ni sustituirlas por veinte helpers nuevos.

**Propuesta:** cargar los checkpoints necesarios y calcular las predicciones; eliminar capas repetidas de firmas, comparaciones y metadatos, conservando las rutinas de figura realmente usadas.

### A4 — Alta: duplicación real de matemáticas PN-IQ en `run_fixed_ridge_sweep.m`

[Ver L120-L173](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/sweep/run_fixed_ridge_sweep.m#L120-L173) y [L248-L299](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/sweep/run_fixed_ridge_sweep.m#L248-L299).

La rutina (**355 líneas**) construye manualmente la misma descomposición en I/Q dos veces: una para identificación y otra para predicción. Es esencialmente la operación ya encapsulada por [`buildPNIQFeatures.m`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/sweep/buildPNIQFeatures.m), que también se emplea en el modelo principal y la inferencia.

**Propuesta:** reutilizar la función científica ya existente en ambos lugares, conservando el orden de columnas `selectedColumns` y la rotación exacta. **Comparar numéricamente los resultados** tras este cambio, pues una divergencia aquí sí alteraría Ridge/NMSE.

### A5 — Alta pero condicionada: infraestructura antigua de algoritmo genético

- [`regressorManager.m`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/pn_gmp_comparison/regressorManager.m): **437 líneas**; conserva propiedades de selección genética, `crossoverrate`, `mutationrate`, evaluaciones y métodos `selection`, `crossover`, `mutation`, `regress`, `prepareForSave`, variantes de bases FV/CVS/MP/DDR y población GMP.
- [`Regressor.m`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/pn_gmp_comparison/Regressor.m): **287 líneas**; contiene generación, ordenación y representación legítimas, pero también métodos de mutación de memoria/orden y cruce propios de genética, con amplios bloques de comentarios explicativos.
- [`GMP_createRegressorManager.m`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/pn_gmp_comparison/GMP_createRegressorManager.m): configura `inittype='GMP'`, llama solo a `initialization`, `removerepeated` y `clearRegressors`; no llama a los métodos genéticos.

**Propuesta:** investigar referencias a métodos de evolución y, si nadie los usa, **retirar del flujo activo el código genético heredado**. No eliminar `initialization`, `removerepeated`, `Regressor`, el orden de regresores ni la construcción de la población GMP sin verificar que siguen saliendo exactamente los **673 regresores** y los mismos términos en el mismo orden. Es una oportunidad de limpieza grande pero con mayor riesgo de alterar la base matemática.

### A6 — Alta: preflight gráfico artificial

[`preflightPaperFigureToolchain.m`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/sweep/preflightPaperFigureToolchain.m) tiene **56 líneas** que comprueban archivos, `exist`, crean una carpeta temporal, generan una figura de prueba, exportan TikZ, comprueban el fichero, registran limpiadores `onCleanup` y destruyen las pruebas.

**Propuesta:** quitar su llamada obligatoria de `run_parameter_sweep`. Cuando se exporte una figura, la propia llamada a `matlab2tikz` dará el error nativo pertinente. Si se desea testear el exportador, mantenerlo como test manual separado, no como dependencia del entrenamiento.

### A7 — Media: validaciones de formatos, contratos y presets en código productivo

- [`countModelFLOPs.m`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/complexity/countModelFLOPs.m): **180 líneas**, múltiples schemas aceptados, validaciones de tablas y convención de FLOPs. Conservar el cómputo científico y la convención, valorar si los distintos formatos son realmente usados.
- [`countModelOperations.m#L122-L140`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/pn_gmp_comparison/countModelOperations.m#L122-L140): lleva incrustado un caso de regresión numérica específico `S==100 && K==17 && F_iq==172` dentro de producción. **Mover esa aserción al test correspondiente**, conservando el cálculo de costes.
- [`plotSweepPaperFigure.m#L139-L171`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/sweep/plotSweepPaperFigure.m#L139-L171): once campos opcionales tratados con `isfield` y default individual. Si solo tiene llamados conocidos, configurar directamente las dos o tres variantes que necesita el proyecto.
- [`resolveMeasurementSampleRate.m`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/sweep/resolveMeasurementSampleRate.m): atiende `fs` y `info_signal.fsovs` y comprueba que coinciden. Si **todas las medidas reales de este proyecto** llevan `fs` y el experimento ya lo ha confirmado, utilizar `measurement.fs` directamente. No eliminar la segunda ruta antes de comprobar los MAT históricos que aún usa la comparación.
- [`buildPhaseNormDataset.m#L103-L128`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/pnnn/buildPhaseNormDataset.m#L103-L128): validación de entradas y órdenes. No es el principal problema, aunque se puede acortar. **No cambiar el orden de las 84 características ni las reglas de fase**.
- [`computeSelectedPointSpectra.m#L7-L24`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/sweep/computeSelectedPointSpectra.m#L7-L24): verifica 3+6 columnas y finitud antes del Welch. Se podría simplificar aceptando que un fallo de tamaño dé el error de MATLAB.
- [`exportPaperFigure.m`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/toolbox/sweep/exportPaperFigure.m): varias comprobaciones de argumentos/rutas antes de `savefig`, `exportgraphics` y `matlab2tikz`. Simplificable, pero su restauración de visibilidad y cierre de figuras sí tiene una finalidad concreta.

### A8 — Media: código de investigación auxiliar con rutas rígidas / documentación desactualizada

- [`run_widely_linear_gmp_probe.m#L66-L84`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/run_widely_linear_gmp_probe.m#L66-L84) carga obligatoriamente `results/parameter_sweep/sweep_d113e389ab78/complexity_sweep.csv`. **Esa carpeta no está en el árbol de `main` auditado** (aparecen `sweep_1f0ac698f55c` y `legacy/sweep_4329deabf42b`). El probe no es autónomo sobre el repositorio publicado; debe apuntar a una fuente explícita, o archivarse como experimento auxiliar.
- [`count_dense_parameters.m`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/count_dense_parameters.m) recompone la población y escribe literalmente `340` en el informe (L83-L84). Si se utiliza de forma habitual, que tome `P` o la configuración; si no, es script exploratorio, no núcleo del DPD.
- [`docs/PLAN_DPD_TRANSACTIONS_2027.md`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/docs/PLAN_DPD_TRANSACTIONS_2027.md) conserva casillas A1–A3 sin marcar y describe las funciones ONLINE/OFFLINE como propuestas, pese a existir y haber sido probadas (L97-L138). Actualizar, no convertir la planificación en lógica ejecutable.
- [`docs/DPD_ILC.md`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/docs/DPD_ILC.md) documenta tanto la nueva interfaz como `testDPD` histórico. Está bien conservar historia si se diferencia claramente de la ruta recomendada.
- [`docs/auditoria_tecnica_repositorio_PNIQ_GMP_PNNN.md`](https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison/blob/main/docs/auditoria_tecnica_repositorio_PNIQ_GMP_PNNN.md) recoge observaciones científicas y del pipeline anteriores; no debe confundirse con una auditoría de simplicidad actual.

### A9 — Media: peso de los resultados versionados

**195/275** archivos versionados están en `results/` y su tamaño agregado en el árbol Git es aproximadamente **53 MB**, sin contar el contenido real LFS. Hay figuras `.fig`, `.png`, `.tikz`, PDFs, CSV y decenas de MAT. Solo cuatro gráficos `.tikz` superan 5,9 MB cada uno. El `.gitignore` ya excluye `/results/dpd_offline/`, pero no evita nuevos artefactos de `results/parameter_sweep/` ni modifica archivos que ya estuvieran versionados.

**Propuesta:** decidir explícitamente qué resultados son artefactos científicos de referencia que deben conservarse y cuáles son regenerables. No borrar retrospectivamente un benchmark histórico ni reescribir Git LFS sin permiso del grupo. Para ejecuciones nuevas, excluir resultados temporales; publicar solo figuras/tablas realmente seleccionadas si esa es la política del proyecto.

## 3. Revisión individual de los 68 archivos MATLAB

Leyenda: **A** intervención prioritaria de simplificación; **B** simplificación localizada; **C** correcto o complejidad científica justificada; **H** histórico/auxiliar; **T** test (el número de aserciones no afecta la ejecución normal).

### 3.1 Raíz y configuración: 10 archivos

| Archivo | Estado | Diagnóstico y acción |
|---|---|---|
| `config/getFairDOMPComparisonConfig.m` | B | 94 líneas de parámetros reales. Mantener presets científicos; revisar nombres/opciones de protocolos antiguos solo al eliminar sus consumidores. |
| `count_dense_parameters.m` | H | Script exploratorio de 84 líneas; imprime P340 fijo. Simplificar entrada o archivar si ya no se usa. |
| `llamadas_DPD.m` | C | Entrada simple para tutor. 46 líneas con espacios, mensajes y ruta de modelo guardado; no sobreingeniería. |
| `main_sweep_and_comparison.m` | C | 33 líneas; orquesta barrido, selección y comparación. La comprobación de override es secundaria. |
| `runDPDOnline.m` | C | 24 líneas físicas (13 sin blancos según informe de Codex). Correcto: load → inferencia → save. |
| `run_parameter_sweep.m` | **A** | 506 líneas; doble identidad SHA, checkpoints hiperfirmados, validaciones del output, transacciones CSV/TXT y preflight obligatorio. Conservar cómputos y checkpoints. |
| `run_selected_comparison.m` | **A** | 525 líneas; cargas/checks repetidos y verificación por hashes, mezclados con figuras útiles. Limpiar control, no borrar gráficas. |
| `run_widely_linear_gmp_probe.m` | B/H | 100 líneas. Experimento científico adicional, pero depende de una ruta de sweep inexistente en `main`. |
| `testDPD.m` | H | 166 líneas del flujo histórico `reuse/recover/train`, no núcleo nuevo. Archivar si se quiere; preservar reproducibilidad. Sigue versionado en GitHub. |
| `trainDPDOffline.m` | C | 59 líneas físicas; flujo limpio y correcto. La complejidad sobrante que sigue usando está *dentro de sus dependencias PNNN*. |

### 3.2 Complejidad, DOMP y métricas: 4 archivos

| Archivo | Estado | Diagnóstico y acción |
|---|---|---|
| `toolbox/complexity/countModelFLOPs.m` | B | 180 líneas; cálculo de FLOPs relevante, pero validaciones exhaustivas de schema/convenio y conversión de formatos a revisar. |
| `toolbox/complexity/getFLOPConvention.m` | C | 26 líneas; tabla de convención FLOP. Mantener. |
| `toolbox/domp/selectDOMPSupport.m` | C | 70 líneas; algoritmo DOMP y least squares. Conservar, no es sobreingeniería. |
| `toolbox/metrics/nmseComplexDb.m` | C | 21 líneas; NMSE y dos condiciones matemáticas de señal válida/energía no nula. No prioritario. |

### 3.3 Generación de bases GMP/PN-IQ: 9 archivos

| Archivo | Estado | Diagnóstico y acción |
|---|---|---|
| `toolbox/pn_gmp_comparison/GMP_createRegressorManager.m` | B | 28 líneas; construye población GMP mediante clase genérica y define parámetros genéticos que el flujo no ejecuta. Candidato a aligerar después de estudiar dependencias. |
| `toolbox/pn_gmp_comparison/Regressor.m` | **A/H** | 287 líneas; la representación y normalización de términos son necesarias; mutaciones/cruces y comentarios de genética parecen legado. Comprobar consumidores antes de limpiar. |
| `toolbox/pn_gmp_comparison/buildGMPRegressorRows.m` | C | 110 líneas de evaluación eficiente por taps periódicos y exponentes; las funciones internas representan subtareas matemáticas claras. Mantener. |
| `toolbox/pn_gmp_comparison/countModelOperations.m` | B | 259 líneas de recuento analítico. Mantener el cálculo; mover el `BaselineRegression` específico de L129-L140 a tests. |
| `toolbox/pn_gmp_comparison/factorizeGMPRegressor.m` | C | 93 líneas de factorización/firmas de términos de la base. Matemática real, mantener. |
| `toolbox/pn_gmp_comparison/modelsConfig.m` | C/H | 30 líneas de constantes históricas para construir diferentes poblaciones; conservar mientras lo requiera `regressorManager.initialization`. |
| `toolbox/pn_gmp_comparison/regressorManager.m` | **A/H** | 437 líneas de clase multiuso con variaciones GMP/MP/DDR/Volterra y genética. El experimento activo solo usa un subconjunto. Reducir con comprobación de población exacta. |
| `toolbox/pn_gmp_comparison/sel_indices.m` | C | 28 líneas; algoritmo de selección de muestras. Mantener. |
| `toolbox/pn_gmp_comparison/selectXYByMapping.m` | B | 29 líneas; admite dos convenciones de mapeo. Si el repositorio solo usa `xy_forward`, se puede sustituir por `x=S.x; y=S.y`, una vez confirmadas dependencias antiguas. |

### 3.4 Red neuronal PNNN: 11 archivos

| Archivo | Estado | Diagnóstico y acción |
|---|---|---|
| `toolbox/pnnn/buildExperimentSignature.m` | **A** | 98 líneas de hashes Java, JSON, bytes y endianness. No aporta matemática PNNN; pertenece al gestor de checkpoints del sweep. |
| `toolbox/pnnn/buildNetworkSignature.m` | **A** | 63 líneas de hash Java de learnables/normalización. Afecta al nuevo OFFLINE indirectamente. |
| `toolbox/pnnn/buildPhaseNormDataset.m` | B/C | 129 líneas; la construcción de features, rotación y taps periódicos es central. Se pueden abreviar 26 líneas de validadores, sin alterar el orden ni `featMode`. |
| `toolbox/pnnn/countPNNNFLOPs.m` | C | 73 líneas de convención de operaciones de la PNNN. Mantener. |
| `toolbox/pnnn/countPNNNParameters.m` | C | 18 líneas de conteo de parámetros; correcta. |
| `toolbox/pnnn/countSparsePNNNFLOPs.m` | C | 47 líneas; conteo sparse y biases, científicamente relevante. |
| `toolbox/pnnn/fit_sparse_pnnn_target.m` | **B/A** | 88 líneas; dos hashes de red densa en L10/L30-L33. Mantener pruning, predicción, resultados y FLOPs. |
| `toolbox/pnnn/predictPhaseNorm.m` | C | 26 líneas de inferencia, directa y necesaria. |
| `toolbox/pnnn/prepare_pnnn_dense_source.m` | B/C | 245 líneas, pero entrena red provisional, selecciona épocas, red final y duración de ajuste: fases científicas reales. Quitar solo firma/validaciones duplicadas, no entrenamientos. |
| `toolbox/pnnn/pruneAndFineTunePNNN.m` | C | 267 líneas de poda global, protección de biases, máscaras de gradiente y Adam por minibatches. Complejidad científica necesaria; no borrar máscaras o controles de integridad de poda. |
| `toolbox/pnnn/summarizeTrainableParameters.m` | C | 48 líneas de máscaras/learnables y recuentos. Mantener, simplificación secundaria. |

### 3.5 Splits y operaciones de sweep: 16 archivos

| Archivo | Estado | Diagnóstico y acción |
|---|---|---|
| `toolbox/splits/buildCommonComparisonSplit.m` | C | 57 líneas de selección reproducible, estratos de amplitud y RNG. Algoritmo real, mantener. |
| `toolbox/sweep/applyDPDModels.m` | B/C | 65 líneas. Inferencia de 3 modelos correcta; contiene compatibilidad opcional con ruta `bundle` y guardado directo. Se puede simplificar la interfaz *solo si* no rompe `testDPD` histórico. |
| `toolbox/sweep/buildPNIQFeatures.m` | C | 36 líneas de features; útil, no eliminar. Reutilizar donde `run_fixed_ridge_sweep` duplica su cálculo. |
| `toolbox/sweep/compactGMPRegressors.m` | C | 10 líneas; compacta la base seleccionada para inferencia. Correcta. |
| `toolbox/sweep/computeSelectedPointSpectra.m` | B/C | 59 líneas; PSD Welch real y contrato de 3+6 columnas; simplificar validaciones si se desea. |
| `toolbox/sweep/exportPaperFigure.m` | B | 61 líneas; guardado FIG/PNG/TikZ con prevalidaciones y restauración de visibilidad. Mantener exportación, reducir defensas. |
| `toolbox/sweep/exportSelectedSpectrumFigures.m` | C | 175 líneas; varios paneles y dos configuraciones de espectros. Largas por el producto gráfico, no por hashes/checks. |
| `toolbox/sweep/fit_complex_gmp_domp.m` | C | 125 líneas de construcción de regresores, DOMP, LS, NMSE, FLOPs y modelo. Mantener. |
| `toolbox/sweep/fit_pniq_gmp.m` | C | 224 líneas de dos proyecciones I/Q, DOMP, ajuste y rotación. Matemática central, no tratar como boilerplate. |
| `toolbox/sweep/getIEEEPaperStyle.m` | C | 22 líneas de estilo visual centralizado. Mantener. |
| `toolbox/sweep/plotSweepPaperFigure.m` | B | 172 líneas de gráficos; `applyDefaults` contiene múltiples `isfield` y variantes configurables. Se puede recortar si solo hay llamadores conocidos. |
| `toolbox/sweep/preflightPaperFigureToolchain.m` | **A** | 56 líneas para crear una figura temporal y testear TikZ; no debe bloquear el entrenamiento numérico. |
| `toolbox/sweep/resolveMeasurementSampleRate.m` | B | 26 líneas comparando `fs` con `info_signal.fsovs`. Sustituir por acceso directo si se confirma un único formato real. |
| `toolbox/sweep/run_fixed_ridge_sweep.m` | **A/B** | 355 líneas de dos familias y varios lambdas. Trabajo numérico real, pero dos bloques PN-IQ duplican `buildPNIQFeatures`; reformar solo esa parte. |
| `toolbox/sweep/run_linear_sweep.m` | C | 52 líneas; coordina las dos familias lineales y soporte nuevo/histórico. Correcta. |
| `toolbox/sweep/selectOperatingPoint.m` | B/C | 218 líneas; selección estabilizada y análisis de sensibilidad son reales; sobran algunas validaciones de tablas/presets, no la lógica de selección. |

### 3.6 Tests: 18 archivos

**No borrar tests solo porque contienen muchos `assert`.** No complican las llamadas normales y sirven para demostrar que una limpieza matemática no cambia resultados. La limpieza de tests es secundaria frente a producción.

| Archivo | Estado | Diagnóstico y acción |
|---|---|---|
| `tests/run_domp_unit_test.m` | T | 90 líneas. Contrasta DOMP con implementación de referencia; valioso. |
| `tests/run_dpd_model_reuse_test.m` | T/H | 102 líneas. Protege reutilización y compatibilidad histórica; conservar mientras exista `testDPD`. |
| `tests/run_dpd_offline_online_test.m` | T | 59 líneas; prueba nueva de OFFLINE/ONLINE. Mantener. |
| `tests/run_dpd_online_fresh_session_test.m` | T | 17 líneas. Verifica ONLINE sin entrenar en sesión nueva. Mantener. |
| `tests/run_experiment_signature_test.m` | T/H | 48 líneas. Quedará obsoleto si se retira SHA-256 y su protocolo; entonces eliminar el test junto con la función. |
| `tests/run_export_paper_figure_roundtrip_test.m` | T | 68 líneas. Test real de FIG/PNG/TikZ; mejor ejecutarlo al preparar el artículo, no antes de entrenar. |
| `tests/run_fixed_lambda_linear_sweep_test.m` | T | 199 líneas. Muchas aserciones sobre protocolo Ridge; proteger resultados si se reutiliza `buildPNIQFeatures`. |
| `tests/run_linear_complexity_sweep_test.m` | T | 168 líneas. Valida soportes, P y FLOPs. Importante para la Transactions. |
| `tests/run_main_sweep_and_comparison_test.m` | T | 73 líneas. Comprueba selección y override sin barrido completo, aunque usa mocks. |
| `tests/run_max_abs_real_parameter_test.m` | T | 154 líneas. Comprueba parámetros normalizados y rango; conservar. |
| `tests/run_measurement_sample_rate_test.m` | T | 31 líneas. Si se adopta `fs` único, ajustar/eliminar ramas antiguas coherentemente. |
| `tests/run_operating_point_selection_test.m` | T | 95 líneas. Test de selección estabilizada; conservar. |
| `tests/run_paper_figure_contract_test.m` | T | 203 líneas. Test de detalles de estilos y formatos; podría reducirse después, sin prioridad. |
| `tests/run_pnnn_shared_dense_sweep_test.m` | T | 171 líneas. Verifica misma fuente densa por presupuesto; conservar especialmente si se elimina hashing. |
| `tests/run_scientific_pipeline_smoke_test.m` | T | 66 líneas. Test rápido de regresores, split, rotación y predicción; conservar. |
| `tests/run_selected_point_spectrum_test.m` | T | 41 líneas. Control del PSD y dimensiones; conservar. |
| `tests/run_sweep_configuration_test.m` | T | 54 líneas. Varias aserciones de constantes del preset; opcional reducir si elimina schemas históricos. |
| `tests/run_sweep_presentation_test.m` | T | 82 líneas. Control de tablas y figura; no prioritario. |

## 4. Resto del repositorio: análisis por tipo, sin fingir lectura de binarios

- `.gitignore`: actualmente ignora `results/dpd_offline/`, no los nuevos resultados de `parameter_sweep`. Decidir política; recordar que archivos previamente trackeados siguen trackeados.
- `.gitattributes`: usa **Git LFS para todos los `.mat`**.
- `.gitmodules`: dependencia de `third_party/matlab2tikz`, no es código propio y su contenido externo no se auditó aquí.
- 108 archivos `.mat`: en el árbol consultado aparecen como punteros LFS de ~131-134 bytes, no sus matrices reales. No permiten deducir si las medidas tienen los campos esperados ni validar redes.
- 195 archivos en `results/`: resultados/figuras, no funciones MATLAB de ejecución. Revisar política de versionado, no suprimir sin distinguir históricos de regenerables.
- Los 13 PDF, 22 FIG, 22 PNG, 22 TikZ, 10 CSV, 2 TXT y TEX son documentación/artefactos, no fuentes `.m`. Se revisaron inventario y rutas, no la integridad interna de todos los contenidos binarios.
- `docs/DPD_ILC.md`, `docs/PLAN_DPD_TRANSACTIONS_2027.md`, auditoría técnica MD/TEX/PDF, reporte PN-IQ PDF y `measurements/README.md`: documentación; la del plan necesita sincronizar estado y los documentos históricos deben etiquetarse como tales. La descripción de derechos de redistribución de medidas en `measurements/README.md` merece respeto al compartir capturas nuevas.

## 5. Plan de limpieza propuesto (sin tocar los experimentos antiguos)

**Fase 1 — Seguridad de trabajo (poco riesgo):** quitar de `run_parameter_sweep` el preflight gráfico obligatorio, sustituir escritura atómica de CSV/TXT por `writetable`/`writelines`, limpiar validaciones de protocolo obsoletas y redundancia del `result`, actualizar documentación, corregir ruta hard-coded del probe. Ejecutar tests de sweep y figuras.

**Fase 2 — Identidades y SHA (riesgo medio):** eliminar `buildExperimentSignature` y `buildNetworkSignature` y sus consumidores **de manera coordinada**, si aceptamos no reanudar automáticamente checkpoints viejos. Preservar `save`/`load`, checkpoint por presupuesto y rutas separadas. Migrar `sweepIdentity` a una configuración simple guardada en el MAT o crear carpetas nuevas para experimentos nuevos. Testear nuevas medidas y comparar inferencia del P340 actual. **No tocar ni reescribir MAT históricos.**

**Fase 3 — Duplicación PN-IQ (riesgo científico controlado):** reemplazar dos bloques manuales de `run_fixed_ridge_sweep` por `buildPNIQFeatures`; comprobar matrices, soportes, coeficientes, NMSE/FLOPs de referencia y figuras antes/después.

**Fase 4 — Legado genético (riesgo científico mayor):** hacer inventario exacto de métodos usados de `Regressor` y `regressorManager`. Tras demostrar que el flujo actual no ejecuta mutación/cruce ni variantes de inicialización, recortar lo heredado. La prueba de aceptación incluye misma población GMP (673), misma ordenación y, si se reentrena una muestra sintética, mismos soportes y predicciones.

**Fase 5 — Acabado:** limpiar presets de opciones que nadie usa y validaciones menores; no convertir los tests en código de producción. No gastar tiempo recortando un `nmseComplexDb` de 21 líneas mientras sigue habiendo un gestor de identidad de cientos.

## 6. Comprobaciones mínimas para aceptar cualquier cambio

1. OFFLINE desde medidas `.mat` de la convención del tutor, para una **captura sintética pequeña**, entrena tres modelos y conserva P.
2. ONLINE carga el `models.mat` real P340 existente y reproduce **idénticamente** las tres señales; no entrena.
3. Nuevo sweep sintético se ejecuta y reanuda sin repetir etapas ya guardadas.
4. La construcción GMP devuelve 673 regresores en el mismo orden; PN-IQ devuelve los mismos índices/features/coeficientes tras la desduplicación.
5. PNNN usa la misma selección de épocas, normalización, entrenamiento denso, máscara y fine-tuning; no confundir quitar hash con cambiar aleatoriedad o método.
6. FLOPs, NMSE, número real de parámetros y lambda de todas las filas permanecen idénticos en fixtures conocidos.
7. Exportación de figuras se prueba **solo cuando se solicitan figuras**, y no bloquea entrenamiento numérico por depender de `matlab2tikz`.

**Límites de verificación:** revisión estática de la rama pública `main`, sin MATLAB ni acceso a archivos locales no subidos. Ninguna propuesta es prueba de equivalencia numérica por sí sola. No se han realizado commits, pushes, modificaciones científicas ni ensayos RF.

## 7. Dictamen

El nuevo OFFLINE/ONLINE ya está bien simplificado a nivel de interfaz, pero **no está completamente libre de la vieja infraestructura**, porque depende de firmas de red basadas en Java a través de `prepare_pnnn_dense_source` y `fit_sparse_pnnn_target`. El mayor volumen de código potencialmente evitable está en el gestor histórico del sweep y en clases GMP heredadas. En cambio, `pruneAndFineTunePNNN`, los ajustes GMP/PN-IQ, las rutinas DOMP, la construcción de regresores y buena parte del código de métricas/figuras reflejan trabajo científico real y no deben reducirse indiscriminadamente.

**Recomendación:** limpiar por etapas pequeñas empezando por A1/A2/A3/A6, no ordenar a Codex «recorta el 70 % de todo el repositorio». Una reducción masiva de las rutinas numéricas podría cambiar de forma silenciosa lo que se está midiendo para el artículo.
