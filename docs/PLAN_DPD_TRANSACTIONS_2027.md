# Hoja de ruta — DPD con GMP-DOMP, PN-IQ-GMP y PNNN

> **Documento de seguimiento del proyecto** · Versión inicial: 11/10/2026  
> **Objetivo a corto plazo:** trabajar con nuevas medidas del tutor mediante un flujo **OFFLINE → ONLINE** sencillo, reproducible y separado.  
> **Objetivo científico:** preparar la comparación **DPD real** para una posible *IEEE Transactions* en enero de 2027.  
> **Estado:** planificación. **No autoriza ejecutar Codex, entrenamientos extensos, transmisiones RF, commits ni pushes.**

## 0. Cómo utilizar este documento

- `[x]` = realizado y respaldado por resultados o código existente.
- `[ ]` = pendiente. No marcarlo como realizado solo porque esté programado o porque pase un test sintético.
- **Prioridad A:** necesario para trabajar con el tutor y sus próximas medidas.
- **Prioridad B:** necesario para defender científicamente los resultados DPD.
- **Prioridad C:** extensión para la Transactions; **no ejecutar todavía**.
- Después de cada sesión, actualizar la tabla de seguimiento del §12 con **fecha, captura, acción, resultado, problemas y siguiente paso**.
- La rama de **modelado del PA** se conserva como referencia histórica; el foco del nuevo artículo será **DPD**.

---

## 1. Qué queremos construir

### Flujo normal (cada vez que el tutor entregue medidas nuevas)

```text
          CAPTURA NUEVA (.mat, x/y)
                      |
                      v
            OFFLINE: ENTRENAR
          +---------------------+
          | GMP-DOMP            |
          | PN-IQ-GMP           |
          | PNNN densa + poda   |
          +---------------------+
                      |
                      v
          GUARDAR 3 MODELOS
      (pesos, coeficientes, configuración,
       normalización y procedencia)
                      |
                      v
           ONLINE: INFERENCIA
            (sin entrenar)
                      |
                      v
       3 SEÑALES DE SALIDA EN .mat
             dpd(k).yvalmod
                      |
                      v
   POSTERIOR MEDIDA RF POR EL AMPLIFICADOR
         (solo en laboratorio autorizado)
```

**Regla esencial:** con cada captura nueva que deba identificarse, empezar un **nuevo entrenamiento completo**. No reutilizar automáticamente modelos, soportes DOMP, pesos, normalizaciones ni coeficientes de otras capturas. Cargar modelos previamente entrenados está permitido **para inferencia** sobre entradas compatibles. La reutilización de cálculos **dentro de una misma captura y mismo experimento** sí será útil, especialmente en el futuro barrido.

**Terminología:** aquí *online* significa **cargar un modelo y generar predicciones sin entrenar**; no implica que el código sea capaz de procesar un flujo RF en tiempo real.

### Convención de señales (evitar el error más peligroso)

- En los `.mat` tipo PNNN, `x` y `y` significan **entrada y objetivo del bloque modelado**, no necesariamente entrada y salida físicas del PA.
- **Imitación del ILC / entrenamiento DPD:** `x = u` (señal de referencia) y `y = x_ILC` (predistorsión calculada por ILC). Entrenamos `u → x_ILC`.
- **Modelado del PA (experimento anterior):** se aprende *entrada efectiva del PA → salida medida del PA*. Comprobar siempre si la entrada procede de `x`, `xCFR` u otra versión procesada; el código de adquisición aplica CFR antes de `measureADRV`.
- La señal de salida de los modelos DPD **no** es la salida medida del PA. Para evaluar linealización hay que transmitirla y medir el PA.
- `mappingMode = "xy_forward"` significa únicamente `x → y` en el archivo seleccionado; no determina por sí solo la función física.

---

## 2. Qué existe y debe preservarse

| Elemento | Estado | Observación |
|---|---|---|
| Tres arquitecturas: Complex GMP-DOMP, PN-IQ-GMP y PNNN | [x] | Ya implementadas en `toolbox/`. |
| Presupuesto común de parámetros reales | [x] | En P=340: 170 coeficientes complejos GMP, 170 características con dos salidas PN-IQ, o 326 pesos + 14 sesgos PNNN. |
| Barrido original **del modelado PA** `20:10:500` | [x] | 49 presupuestos; **no es todavía el barrido de DPD**. |
| Selección original de P=340 | [x] | Criterio de estabilización conjunta aplicado al modelado PA; **no asumir que es óptimo para DPD**. |
| Primer experimento de imitación ILC con P=340 | [x] | NMSE offline y modelos recuperados. |
| Carga e inferencia programática con `applyDPDModels` | [x] | Trabaja con un `bundle` y vector de entrada; falta entrada cómoda por archivo/flujo de usuario. |
| Recuperación histórica `testDPD(340,"recover")` | [x] | Es recuperación de **la misma captura**, no entrenamiento desde cero de captura nueva. |
| Selección libre de medidas para OFFLINE | [ ] | `testDPD.m` fija un `.mat` concreto. |
| OFFLINE/ONLINE separados y fáciles de ejecutar | [ ] | Falta completar las dos entradas operativas. |
| Barrido específicamente DPD | [ ] | Futuro; no ejecutarlo ahora. |
| Ensayos RF con los tres DPD | [ ] | No inferirlos a partir del NMSE offline. |

**Resultados de referencia para no mezclarlos:**

| P=340 | Modelado PA, NMSE señal completa* | Imitación ILC, NMSE offline* |
|---|---:|---:|
| Complex GMP-DOMP | −38,860 dB | −38,456 dB |
| PN-IQ-GMP | −42,740 dB | −38,931 dB |
| PNNN | −39,243 dB | −35,630 dB |

\* Las dos columnas corresponden a **objetivos diferentes**. Además, la evaluación denominada *señal completa* contiene el subconjunto usado para entrenar: **no es un test independiente**.

**Protección del trabajo anterior:** no sobrescribir `.mat` históricos ni cambiar sus resultados de referencia. Guardar cada entrenamiento nuevo en un directorio específico de captura/configuración/ejecución; conservar la procedencia y, cuando corresponda, hashes de medidas/modelos.

---

## 3. FASE A — Preparar OFFLINE/ONLINE para el laboratorio [PRIORIDAD A]

### A1. Selección explícita de medidas

- [ ] Permitir que el usuario proporcione **ruta absoluta o relativa** a un `.mat` de entrenamiento.
- [ ] Mantener el contrato del repo PNNN: campos `x` y `y`; metadatos opcionales `fs`, `fsmed`, `description` e `info_signal` cuando existan.
- [ ] Validar existencia, dimensiones, mismo número de muestras, tipo numérico compatible y ausencia de NaN/Inf.
- [ ] Exponer y guardar la **interpretación de x/y** y la identidad de la captura; no deducir el bloque físico de `xy_forward`.
- [ ] Conservar frecuencia de muestreo y escalas físicas conocidas; **no inventarlas** si faltan.
- [ ] Verificar la convención exacta con el tutor si las nuevas capturas contienen `meas_out` en lugar de `x/y` directas.

**Aceptación A1:** se pueden seleccionar **dos archivos distintos** sin modificar el código fuente ni editar rutas incrustadas.

### A2. OFFLINE: entrenamiento completo de los tres modelos

- [ ] Crear una entrada simple, por ejemplo `trainDPDOffline(measurementFile,P)`; el nombre es **propuesta, no API existente**.
- [ ] Para la captura seleccionada, ejecutar la selección DOMP **nueva** de cada modelo lineal y estimar sus coeficientes.
- [ ] Inicializar y entrenar **desde cero** una PNNN densa; podar al presupuesto P y realizar el reajuste posterior.
- [ ] Impedir *warm start* accidental, carga silenciosa de otro experimento o uso de normalización ajena. **El antiguo PNNN tiene opciones de warm start: no copiarlas como comportamiento obligatorio.**
- [ ] Mantener el protocolo científico actual del repositorio, arquitectura de PNNN, semillas y presupuesto real P; **no sustituirlos** por la configuración por defecto del PNNN original.
- [ ] Entrenar **un solo presupuesto P** por ejecución normal. **No** arrancar `20:10:500`.
- [ ] Guardar individualmente coeficientes/regresores/soportes GMP; soporte I/Q y coeficientes PN-IQ; red **final podada** y normalización PNNN.
- [ ] Guardar bundle/deploy autocontenido, resumen de NMSE/FLOPs/parámetros, configuración, firma de captura y tiempo de entrenamiento.
- [ ] No sobrescribir sesiones previas. Una recuperación tras fallo debe ser **explícita**, limitada a la misma identidad de experimento y distinguible de entrenar desde cero.

**Aceptación A2:** la captura A produce tres modelos nuevos y un paquete reutilizable; la captura B genera **otro conjunto independiente**, aunque se use el mismo presupuesto P.

### A3. ONLINE: inferencia sin entrenamiento

- [ ] Crear una entrada simple, por ejemplo `runDPDOnline(modelFile,inputFile)`; nombre **propuesto, no API existente**.
- [ ] Seleccionar explícitamente **qué paquete de modelos** se carga y **qué archivo de señal** se procesa; leer `x` (o alias documentado) con la lógica compatible de PNNN.
- [ ] Reutilizar `applyDPDModels` para calcular las tres salidas; **no** ejecutar DOMP, `trainnet`, poda o *fine-tuning*.
- [ ] Reutilizar las estadísticas de normalización del entrenamiento; aplicar a la nueva entrada solo el tratamiento previsto (p. ej., centrado de su propia media cuando así quedó definido).
- [ ] Comprobar compatibilidad de frecuencia de muestreo, amplitud y memoria periódica; rechazar o advertir cuando no haya información suficiente.
- [ ] Guardar `dpd(k).yvalmod`, `dpd(k).modeltype` (3 modelos), más metadatos de captura/modelo/configuración.
- [ ] La salida debe poder generarse aunque no haya referencia `y` para evaluar NMSE.
- [ ] Si existe objetivo medido, calcular NMSE **opcional** y sin mezclarlo con medidas RF posteriores.
- [ ] No sobreescribir señales de otras ejecuciones; probar con MATLAB recién iniciado.

**Aceptación A3:** cerrar MATLAB, abrir una sesión nueva, seleccionar el deploy A y una señal compatible y obtener tres vectores complejos **sin entrenamiento**.

### A4. Pruebas antes de ir al laboratorio

- [ ] Test sintético ligero: dos capturas diferentes → dos identidades de entrenamiento distintas.
- [ ] Comprobar recuento **exacto** de parámetros P en cada modelo.
- [ ] Comprobar que guardar/abrir el deploy no cambia las predicciones.
- [ ] Confirmar mediante tests o instrumentación que ONLINE **nunca** llama al entrenamiento.
- [ ] Verificar la longitud de señal, orden de modelos, `yvalmod`, formato `.mat`, no sobrescritura y registro de procedencia.
- [ ] Ejecutar un OFFLINE real **solo con autorización expresa** (es costoso).
- [ ] Generar ONLINE con la misma captura de entrenamiento y, si se dispone, con una segunda señal compatible.
- [ ] Revisar compatibilidad de las salidas con la cadena de medida: media/DC, escalado, límites DAC, CFR/PAPR, longitud y sincronización.

**Hito A — Listo para laboratorio:** el tutor puede entregar `medida.mat`; se completan **OFFLINE → ONLINE** sin tocar el código ni volver a entrenar durante ONLINE.

---

## 4. FASE B — Validación científica DPD rigurosa [PRIORIDAD B]

### B1. Separar entrenamiento y evaluación

> Crítica del Reviewer 2: usar el ~10 % para identificar y reportar NMSE sobre el 100 % mezcla muestras vistas y no vistas.

- [ ] Mantener la partición **histórica** como referencia reproducible, sin modificar resultados pasados.
- [ ] Definir una evaluación principal sobre el **~90 % no empleado para ajustar** la versión final de cada modelo; documentar exactamente índices y posible contacto temporal con muestras de entrenamiento por efectos de memoria.
- [ ] Conservar un conjunto **interno** de validación si se necesitan seleccionar épocas/ajustes, sin contaminar el conjunto externo de prueba.
- [ ] Para las tres arquitecturas, utilizar las **mismas muestras de prueba**, misma definición de NMSE y escalado.
- [ ] Distinguir en tablas/figuras: `NMSE_identificación`, `NMSE_test_90%`, `NMSE_señal_completa` (solo diagnóstico) y, cuando exista, `NMSE_RF_salida_PA`.
- [ ] Considerar pruebas con **otra captura independiente** para estudiar generalización a otros estados del PA.

**Aceptación B1:** nunca llamar *test independiente* a la señal completa de la que proceden las muestras de entrenamiento.

### B2. Medir la calidad real de los predistorsionadores

- [ ] Obtener salidas de los tres DPD entrenados para la **misma señal deseada y estado del PA**.
- [ ] Medir sobre el mismo montaje, potencia y procedimiento: referencia sin DPD, ILC de referencia (si procede), GMP-DOMP, PN-IQ-GMP y PNNN.
- [ ] Registrar **entrada realmente enviada al PA** (posible CFR), salida PA capturada, sincronización, ganancias, amplitudes y límites del equipo.
- [ ] Calcular **NMSE de salida PA** respecto de la referencia adecuada; **ACPR/ACLR**, espectros y **EVM** cuando sea aplicable.
- [ ] Guardar señales y metadatos de RF sin sustituir las inferencias originales.
- [ ] Garantizar comparabilidad (misma waveform, frecuencia de muestreo, potencia, configuración del amplificador y procedimiento de adquisición).

**Aceptación B2:** tabla experimental de linealización; **no** presentar NMSE de imitación ILC como NMSE físico del PA.

### B3. Fiabilidad numérica y física

- [ ] Investigar la característica de entrada PNNN teóricamente nula que tiene `sigmaX` casi cero (~8×10⁻¹⁹ en el checkpoint histórico). No modificar el baseline silenciosamente: hacer variante controlada con nueva validación si procede.
- [ ] Aclarar cómo se normalizan/acondicionan las columnas de cada familia y qué magnitudes son comparables entre métodos.
- [ ] Documentar coste de generación de características, retardos, operaciones especiales y condiciones para aprovechar la poda en ejecución real.
- [ ] Analizar físicamente los resultados cuando exista evidencia (memoria, AM/PM, amplitud y efectos de potencia), evitando explicaciones no demostradas.

---

## 5. FASE C — Barrido DPD de complejidad para la Transactions [PRIORIDAD C; NO EJECUTAR AHORA]

**Objetivo:** repetir el análisis `P = 20:10:500` pero para la relación **`u → x_ILC`** y, donde se mida el PA, para la eficacia de linealización.

- [ ] Diseñar una entrada `runDPDParameterSweep(measurementFile,parameterGrid)` (**API propuesta**), que **reutilice funciones OFFLINE**, no una implementación independiente.
- [ ] Una ruta DOMP por arquitectura lineal y captura; usar sus prefijos para varios presupuestos de **la misma captura**.
- [ ] Una red densa PNNN por captura y configuración; generar podas independientes y *fine-tuning* por presupuesto sin entrenar la densa desde cero 49 veces.
- [ ] Guardar checkpoint por presupuesto con identidad y controles estrictos; reanudar explícitamente si se interrumpe.
- [ ] Guardar **modelos inferibles**, métricas y metadatos por punto; no limitarse a predicciones.
- [ ] Para cada presupuesto, calcular parámetros reales activos, FLOPs, NMSE offline de identificación/test y coste numérico adicional pertinente.
- [ ] Exportar curvas: NMSE vs parámetros, NMSE vs FLOPs, FLOPs vs parámetros, sparsity vs calidad.
- [ ] Seleccionar presupuestos representativos para ensayos RF y añadir **ACLR/ACPR/EVM vs parámetros/FLOPs** (si hay suficientes mediciones).
- [ ] Seleccionar **un nuevo punto operativo específico para DPD**; no trasladar automáticamente el P=340 escogido en modelado PA.

**Regla de rendimiento:** el barrido debe poder reutilizar cálculos **dentro de una captura**, sin reutilizar modelos de otra campaña.

---

## 6. FASE D — Comparaciones de complejidad y calidad [PRIORIDAD C]

| Comparación | Pregunta | Condición para hacerlo bien |
|---|---|---|
| **Igual número de parámetros reales activos** | ¿Quién aprovecha mejor P coeficientes/pesos? | Misma definición de parámetro y mismo conjunto de datos. |
| **Igual presupuesto de FLOPs** | ¿Quién funciona mejor con un coste aritmético parecido? | Elegir puntos bajo un techo común de FLOPs o interpolar curvas justificándolo; los FLOPs reales pueden ser discretos. |
| **Igual calidad objetivo** | ¿Quién alcanza cierto NMSE / ACLR / EVM con menos parámetros o FLOPs? | Calidad evaluada con la **misma métrica y el mismo conjunto de prueba**; RF solo cuando exista medida. |
| **Coste de implementación real** | ¿Quién consume menos recursos, latencia o energía? | Requiere prototipo/benchmark específico; los FLOPs son solo un indicador. |

- [ ] Conservar una **convención explícita de conteo** (sumas/multiplicaciones reales y complejas; activaciones, divisiones y generación de características reportadas de forma transparente).
- [ ] Distinguir los FLOPs **ideales de inferencia sparse** de PNNN del tiempo real de MATLAB y del hardware.
- [ ] Mantener tablas que permitan comparar por parámetros, FLOPs y calidad sin reentrenar para generar cada figura.
- [ ] Cuando las arquitecturas tengan diccionarios candidatos diferentes, explicitar que **igual P no implica idéntico espacio de representación**.

---

## 7. FASE E — Estudio de ablación de PN-IQ-GMP [PRIORIDAD C]

**Nombre del experimento:** **estudio de ablación (*ablation study*)**. Objetivo: separar la contribución de la **normalización de fase (PN)**, de la representación/ajuste **I/Q** y de su posible interacción.

| Variante conceptual | PN | Representación I/Q independiente | Estado |
|---|:---:|:---:|---|
| GMP de referencia | No | No | [x] Disponible como referencia. |
| PN-GMP | Sí | No | [ ] Definir qué cambio representa realmente. |
| IQ-GMP | No | Sí | [ ] Definir selección y restricciones de coeficientes. |
| PN-IQ-GMP | Sí | Sí | [x] Modelo existente; falta ablación completa. |

**Advertencias científicas clave:**

- Rotar y desrotar la fase de un **GMP complejo lineal** puede ser una transformación **matemáticamente equivalente** al GMP original. No atribuir a PN una mejora que no pueda existir bajo la definición propuesta.
- PN-IQ modifica también el **diccionario disponible**, la **selección del soporte** y la **libertad de los coeficientes**. Hacen falta controles adicionales para aislar cada causa, no solo cuatro nombres de modelos.
- Incluir, cuando esté bien definido, un control por **libertad de coeficientes** y otro por **tamaño/composición del diccionario**, para atender la duda del Reviewer 2 sobre *degrees of freedom*.
- El experimento auxiliar **widely-linear GMP** ya existe para el modelado PA; valorar su reproducción/adaptación a DPD, sin presentarlo como sustituto de toda la ablación.
- Comparar variantes con el **mismo P** y, en una fase posterior, por **coste parecido**. Usar las mismas medidas, splits y criterios.
- [ ] Analizar si hay una **interacción PN × IQ**, es decir, si la combinación aporta más que la suma de sus efectos aislados. Evitar afirmar de antemano que existe sin datos.

**Aceptación E:** una tabla y curvas que permitan explicar, con controles matemáticamente bien definidos, **de dónde procede la ganancia de PN-IQ-GMP**.

---

## 8. FASE F — Revisiones científicas incorporadas al plan

Estas observaciones provienen de los **Reviewer 1 y 2** compartidos para el trabajo anterior. Se consideran requisitos de análisis para la Transactions, no órdenes de modificar sin validar.

| Origen | Observación | Acción | Prioridad | Estado |
|---|---|---|---|---|
| R1 | ¿Qué regularización utiliza la función de coste de PNNN? En la tabla falta valor. | Documentar **explícitamente si no existe regularización L1/L2**; si procede, estudiar una variante regularizada y reportar el coeficiente. No confundir poda con penalización L1/L2. | B/C | [ ] |
| R1 | Añadir ACPR en la tabla. | Medir y publicar **ACPR/ACLR del PA** con las señales DPD; no inventarlo a partir de NMSE offline. | B | [ ] |
| R1 | ¿Mejora la distorsión fuera de banda al aumentar la complejidad PNNN? | Barrido de complejidad DPD y espectros/medidas RF en puntos seleccionados. | C | [ ] |
| R2 | Se identifica con 10 % de ~49k muestras y se valida con el 100 %. | Reservar y reportar el **90 % no usado para entrenamiento** como evaluación externa, además de conservar el histórico por compatibilidad. | B | [ ] |
| R2 | ¿La mejora procede de PN o de más grados de libertad? | Ablación factorial PN/IQ + controles de libertad de coeficientes/diccionario. | C | [ ] |

**Nota de alcance:** comentarios de otras revisiones sobre POMP/Volterra o formulaciones diferentes deben inventariarse **por separado** antes de incorporarlos; no asumir que se refieren a este mismo manuscrito.

---

## 9. FASE G — Redacción y figuras de la Transactions [PRIORIDAD C]

- [ ] Establecer la **pregunta científica central** de la Transactions: eficacia DPD real y relación precisión/complejidad de GMP-DOMP, PN-IQ-GMP y PNNN.
- [ ] Separar claramente resultados de **modelado PA**, **imitación offline del ILC** y **linealización experimental PA**.
- [ ] Explicar el origen del conjunto de medidas: ILC y PA, variables `u`, `x`, `xCFR`, `y`; ruta de ganancia, sincronización y preprocesado.
- [ ] Mostrar comparación justa en parámetros, FLOPs y nivel de calidad, sin confundir FLOPs ideales con implementación real.
- [ ] Incluir ablación PN/IQ y controles de grados de libertad.
- [ ] Incluir pruebas de generalización con evaluación no entrenada y, si se puede, captura independiente.
- [ ] Incluir ACPR/ACLR, EVM y espectros medidos cuando estén disponibles.
- [ ] Presentar incertidumbres y límites: efectos de memoria, sensibilidad a potencia, normalización y recursos de hardware.
- [ ] Mantener trazabilidad tabla/figura → CSV/MAT/captura/modelo/configuración/script.

**Fecha orientativa de trabajo:** enero de 2027. Confirmar convocatoria/plazo real antes de fijar calendario irreversible.

---

## 10. Orden de ejecución recomendado (sin lanzar nada ahora)

1. **Revisar contrato de medidas con el tutor:** archivo ejemplo, x/y, frecuencia, potencia y propósito DPD; confirmar tratamiento de `x` y `xCFR` para comparaciones físicas.
2. **Implementar A1–A3:** selección de medidas, OFFLINE completo para P fijo, deploy y ONLINE por archivo.
3. **Pasar A4:** pruebas sintéticas y una corrida real autorizada; verificar que ONLINE no entrena.
4. **Comprobar uso en laboratorio (B2):** generar los tres DPD, transmitir y medir con controles RF apropiados.
5. **Cerrar B1/B3:** partición independiente, tratamiento numérico y documentación de complejidad.
6. **Después, sin afectar el flujo normal, desarrollar C/D/E:** barrido DPD, FLOPs/igual calidad y ablaciones.
7. **Preparar G:** tablas/figuras reproducibles, interpretación y texto del artículo.

**No hay que implementar todo a la vez.** Cada etapa exige pasar su aceptación antes de darla por cerrada.

---

## 11. Decisiones abiertas — resolver antes de codificar o transmitir

| ID | Pregunta | Situación / criterio |
|---|---|---|
| D1 | ¿El `.mat` nuevo siempre contiene `x`, `y`? | El usuario indica que sí; validar con la primera captura real. |
| D2 | ¿Cómo se selecciona la señal para ONLINE? | Permitir archivo de entrada independiente del de entrenamiento; prioridad a `x` según `xy_forward`. |
| D3 | ¿Puede entrenarse de nuevo la **misma** captura a petición del usuario? | Sí, crear **otra ejecución sin sobrescribir**; la recuperación no debe sustituir al entrenamiento explícito. |
| D4 | ¿Qué datos usa el NMSE principal? | Plan: identificación para ajuste, test no visto para métrica principal; conservar resultado histórico en campo distinto. |
| D5 | ¿Cuál es la entrada física correcta del modelo PA? | Confirmar `x` frente a `xCFR` y el preprocesado de `_forward_xy.mat`; no inferirlo del nombre. |
| D6 | ¿Cuál es el criterio para reentrenar ante cambios del PA? | Por defecto cada nueva captura de identificación → entrenamiento nuevo. Inferencia con modelo antiguo solo como prueba explícita de generalización. |
| D7 | ¿Cómo comparar con presupuesto igualado de FLOPs? | Definir umbrales/discretización y qué operaciones computar antes de publicar resultados. |
| D8 | ¿Qué cuatro variantes exactas constituyen la ablación? | Requieren especificación matemática que evite variantes equivalentes artificiales. |
| D9 | ¿Qué metadatos del banco RF son imprescindibles? | Confirmar ganancia, potencia, fs, alineación, CFR/DAC y métricas con el tutor. |
| D10 | ¿Cómo resolver el canal de desviación casi nula PNNN? | Investigar fuera de baseline; cambio solo con nueva versión experimental. |

---

## 12. Registro de avance (actualizar tras cada sesión)

| Fecha | Captura / versión | Paso | Estado | Evidencia / archivo | Próxima acción |
|---|---|---|---|---|---|
| 2026-10-09 | `experiment20260429T134032_xy.mat` · P340 | Recuperación de modelos DPD históricos | Completada **según validación MATLAB reportada** | `docs/DPD_ILC.md`; archivos locales `models_P340.mat`, `validation_20261009_001045_626.mat` | Diseñar flujo nuevo sin romper histórico |
| 2026-10-11 | Planificación | Definición OFFLINE/ONLINE por archivo y roadmap Transactions | Documentada | Este `.md` | Empezar por A1 cuando se autorice |
| 2026-10-11 | Dos capturas sintéticas, 512 muestras, P20; baseline P340 preservado | Implementación A1–A3 y pruebas ligeras A4 | Código implementado; validación sintética y MATLAB reiniciado superadas. No se da por cerrado el hito de laboratorio | `trainDPDOffline.m`, `runDPDOnline.m`, `docs/DPD_ILC.md`, tests `run_dpd_offline_online_test` / `run_dpd_online_fresh_session_test`; 6 DOMP + 6 trainnet OFFLINE, 0 entrenamiento ONLINE, predicciones idénticas tras recarga | Autorizar primer OFFLINE real; confirmar semántica, fs, amplitud, DC y compatibilidad con el tutor. Sin RF ni barrido DPD |
| — | — | — | Pendiente | — | — |
| 2026-10-11 | Dos tests sintéticos P20 y P340 real ya guardado | Simplificación OFFLINE/ONLINE solicitada | Código de producción 247 → 58 líneas; tres helpers eliminados. Tests y recarga pasan; P340 reproduce exactamente su ONLINE anterior (error relativo 0) | `trainDPDOffline.m`, `runDPDOnline.m`, `docs/DPD_ILC.md`; sin Java de archivos, hashes de medidas ni gestión UUID propia | Mantener formato conocido x/y/fs y protocolo científico; sin reentrenar P340 ni transmitir RF |

### Plantilla breve para próximas sesiones

```text
Fecha:
Archivo de medidas (ruta/nombre):
Objetivo: OFFLINE / ONLINE / RF / SWEEP / ABLACIÓN / OTRO
Modelo(s) y presupuesto(s):
Qué se ejecutó:
Resultados (NMSE identificación/test, FLOPs, ACPR/ACLR, EVM si procede):
Dónde se guardaron los modelos y señales:
Problemas o decisiones tomadas:
Qué toca después:
```

---

## 13. Referencias operativas

- Repositorio actual: `sermunagu/PNNN_PNGMP_DOMP_fair_comparison`.
- Repositorio del flujo original: `sermunagu/PNNN` (`train_PNNN_offline.m`, `run_PNNN_online_from_xy.m`, `config/getPNNNConfig.m`).
- Lógica DPD actual: `testDPD.m`, `toolbox/sweep/applyDPDModels.m`, `docs/DPD_ILC.md`.
- Barrido anterior: `run_parameter_sweep.m`, `run_selected_comparison.m`, `toolbox/sweep/selectOperatingPoint.m`.
- Modelos y datos: `measurements/experiment20260429T134032_xy.mat` para imitación ILC; `_forward_xy.mat` para comparación anterior del PA.
- Código adquisición ILC: `sermunagu/ILC_PAWR2027_v20260428`, donde aparecen las señales `u`, `x`, `xCFR`, `y` y `meas_out`.

**Regla de conservación:** incorporar este documento a `docs/` cuando se decida actualizar el repositorio. Crear el `.md` por sí solo **no modifica el GitHub ni implica ejecutar los experimentos**.
