# Diagnóstico numérico de la normalización de fase PNNN

> Organización actual (14/10/2026): ejecutores PA en `PA/` y DPD en `DPD/`;
> `runPNIQAblation` y `runPNNNPhaseExperiment` están en `DPD/experiments/`.
> Antes de usar ejemplos: `addpath('PA','DPD','PA/experiments','DPD/experiments')`.
> Las narraciones históricas conservan sus nombres; el [mapeo actual](GUIA_REPOSITORIO_DPD.md#entradas-reorganizadas-14102026)
> indica sus sustitutos, sin cambiar medidas, modelos ni resultados.

Fecha: 2026-10-12. MATLAB R2026a. DPD offline por imitación del ILC, sin transmisión RF.

> Actualización 13/10/2026: el entrenamiento experimental P310 de cero
> estructural sí se completó en
> `results/pnnn_phase_experiment/tp1bfa4378_d472_497b_87d7_b5c638c34392/`.
> NMSE de evaluación: histórica −35.482326 dB; experimental −35.525255 dB
> (mejora 0.042929 dB). La declaración de ejecución pendiente más abajo
> pertenece al diagnóstico inicial. No se sustituyó el baseline ni se validó RF.

## Resultado principal

Existe una amplificación real del redondeo en la característica 15: un residuo de desviación típica `8.06e-19` se convierte, mediante z-score, en una entrada de desviación típica `1.00546`. Sin embargo, **no afecta directamente a las predicciones del P310 seleccionado**, porque sus 12 conexiones hacia `fc1` están podadas. P310 mantiene la equivariancia dentro de la precisión numérica de la inferencia.

Sí hay un efecto causal medible en la red densa guardada y, mucho menor, en P500. Esto justifica preparar una variante experimental de cero estructural, no sustituir el baseline ni atribuirle automáticamente la diferencia de NMSE respecto a los modelos lineales.

No se ha ejecutado el entrenamiento real corregido P310 ni un nuevo barrido. Se han realizado diagnóstico con modelos existentes, intervención de sensibilidad con pesos congelados y pruebas sintéticas breves de la implementación experimental.

## A. Diagnóstico e implementación auditada

### Datos y protocolo real

- Repositorio: `C:\Sergi\Investigacion\Códigos\PNNN_PNGMP_DOMP_fair_comparison`.
- Barrido: `results/dpd_sweep/tp98e6f7ac_f3cd_44c7_9ce3_fd8619319748/`.
- Captura: `measurements/experiment20260429T134032_xy.mat`.
- Señales: entrada `x = u`, objetivo `y = x_ILC`; 491520 muestras complejas.
- Identificación: 49151 muestras; entrenamiento interno: 41778. El resto de identificación se utiliza para seleccionar épocas.
- Se conserva la selección histórica de un segmento que contiene el máximo del objetivo. La eliminación de DC usa la captura completa y los taps usan extensión periódica.
- El complemento de identificación pertenece a la misma captura, y además sus curvas se han usado para seleccionar P310. **No es un test independiente** ni una campaña RF.

`trainDPDOffline` y `runDPDParameterSweep` llaman a `prepare_pnnn_dense_source`. Esta función construye las características mediante `buildPhaseNormDataset`, obtiene estadísticas del entrenamiento interno y selecciona la época densa con validación interna. Selecciona también la duración de fine-tuning con el presupuesto de referencia P200. Después descarta la red provisional, reinicializa con la misma semilla y recalcula las estadísticas sobre toda la identificación para el ajuste denso definitivo.

`fit_sparse_pnnn_target` parte de esa fuente densa. `pruneAndFineTunePNNN` mantiene la poda por magnitud, protege los 14 biases, fija las máscaras y vuelve a aplicar los ceros durante el fine-tuning. No se ha modificado ninguna de esas reglas, la semilla 42 ni su offset de fine-tuning. El baseline real seleccionó 1228 épocas densas y 164 de fine-tuning.

La arquitectura sigue siendo `84 -> fc1(12) -> sigmoid -> fcOut(2)`: 1032 pesos y 14 biases en la red densa. `summarizeTrainableParameters` conserva el recuento por máscaras para redes podadas. P310 tiene 296 pesos activos y 14 biases, con 780 FLOPs/muestra según el contador histórico.

ONLINE llama a `applyDPDModels`: elimina DC como el entrenamiento, construye taps periódicos, usa `model.featureConfig` y las estadísticas guardadas, y llama a `predictPhaseNorm`. Esta última desnormaliza las dos salidas y deshace la rotación. No estima normalizaciones nuevas ni entrena.

### Origen de la característica 15

Para `M=13`, modo `full`, el orden es:

| Posiciones | Características |
|---|---|
| 1–14 | Parte real de los 14 taps rotados |
| 15–28 | Parte imaginaria de los 14 taps rotados |
| 29–42 | Envolventes de orden 1 |
| 43–56 | Envolventes de orden 3 |
| 57–70 | Envolventes de orden 5 |
| 71–84 | Envolventes de orden 7 |

Para el tap actual no nulo, `rotation = conj(x)/abs(x)` y, en aritmética exacta,

\[
\operatorname{Im}\{(\overline{x}/|x|)x\}=\operatorname{Im}\{|x|\}=0.
\]

La división y el producto complejos se redondean. El residuo imaginario es del orden de la precisión de máquina por la amplitud de entrada. La protección histórica sustituye por 1 únicamente las desviaciones **exactamente** cero; no reconoce este residuo. Por ello lo divide por una sigma del mismo orden y obtiene una variable de orden unidad que carece de significado físico.

Estadísticas P310, señal completa sin rotación adicional:

| Magnitud | Valor |
|---|---:|
| Media guardada en identificación | `1.567958675115e-21` |
| Sigma guardada en identificación | `8.019388004715e-19` |
| Media original, señal completa | `1.116216252658e-21` |
| Desviación típica original, señal completa | `8.063174048349e-19` |
| Mínimo / máximo originales | `-6.938893903907e-18` / `+6.938893903907e-18` |
| Media después del z-score guardado | `-0.000563313` |
| Desviación típica después del z-score | `1.005460023` |
| Mínimo / máximo después del z-score | `-8.654602893` / `8.650692473` |
| Sigma guardada / RMS de entrada | `3.192575868e-17` |

Estos valores normalizados son artificialmente grandes respecto al cero matemático, aunque no sean los mayores de toda la matriz de características. En particular, el orden 7 alcanza valores normalizados de 19.60 por motivos legítimos.

### Otras escalas: no aplicar un umbral absoluto general

| Grupo | Rango de sigma guardada | Máximo absoluto normalizado |
|---|---:|---:|
| Real de taps, 1–14 | `0.0116430` – `0.0177723` | 4.309 |
| Imaginaria de memoria, 16–28 | `0.00618864` – `0.0178324` | 4.509 |
| Envolvente orden 1 | `0.0116430` – `0.0116458` | 3.799 |
| Envolvente orden 3 | `3.14578e-5` – `3.14714e-5` | 8.666 |
| Envolvente orden 5 | `8.88986e-8` – `8.89339e-8` | 14.223 |
| Envolvente orden 7 | `2.88815e-10` – `2.88890e-10` | 19.601 |

No se identifica otro canal casi nulo por cancelación estructural. Las potencias altas son pequeñas porque la amplitud de entrada es menor que uno, pero tienen dispersión física: en orden 7, `std(feature)/RMS(feature) >= 0.9736`. Comparar directamente sus sigmas con un umbral absoluto común mezclaría escalas y potencias distintas. También existen redundancias teóricas entre características, como real del tap actual y su envolvente de orden 1; no se han eliminado ni modificado.

## B. Evidencia con modelos reales guardados

### Rotaciones y definición del error

Se aplicó `phase*u` al ONLINE real, incluida la eliminación de DC y los taps periódicos. Se reutilizaron exactamente las estadísticas y modelos guardados. Se calculó:

\[
E(\theta)=\frac{\|f(e^{j\theta}u)-e^{j\theta}f(u)\|_2}{\|f(u)\|_2}.
\]

Las normas se calcularon en `double`, convirtiendo cada señal por separado. No se concatenaron las predicciones lineales con la PNNN en una matriz que pudiera reducir su precisión a `single`. Para NMSE se rotó también el objetivo centrado. Cada intervención tiene su propia referencia de 0° para E; el cambio entre baseline e intervención se mide aparte.

| Rotación | Cambio RMS del canal 15 normalizado | E GMP-DOMP P310 | E PN-IQ-GMP P310 | E PNNN P310 |
|---:|---:|---:|---:|---:|
| 0° | 0 | 0 | 0 | 0 |
| 30° | 1.42390 | `1.53895e-3` | `1.80981e-3` | `3.57395e-8` |
| 60° | 1.42173 | `2.66860e-3` | `3.14264e-3` | `3.57459e-8` |
| 90° | 0.96509 | `3.09042e-3` | `3.65358e-3` | `6.16058e-17` |
| 135° | 1.42289 | `2.23061e-3` | `2.69117e-3` | `3.57318e-8` |
| 180° | 1.27179 | `5.75538e-4` | `9.66362e-4` | `1.28426e-16` |
| 270° | 1.35965 | `3.09025e-3` | `3.65082e-3` | `1.88876e-16` |

La matriz completa de entradas normalizadas cambia un 15.55 % en norma Frobenius a 30°, dominada por este canal. Eso **no implica** que cambie la salida: en P310 las 12 conexiones del canal 15 son exactamente cero y sus máscaras son falsas. Su contribución a las 12 preactivaciones es cero; las contribuciones RMS de las otras entradas oscilan entre 0.41264 y 1.5816.

Los controles lineales no son controles de equivariancia exacta: contienen los regresores auxiliares `conj(x)` y `abs(x)`, con grados de fase -1 y 0, en lugar de 1. En los modelos compactos P310 aparecen en posiciones 36/126 (GMP) y 42/80 (PN-IQ). Bajo rotación esos términos no siguen la misma ley que `x*abs(x)^p`. Por tanto sus errores de orden `1e-3` no deben atribuirse al canal 15 ni interpretarse como un fallo numérico equivalente. No se ha cambiado su diccionario, soporte o ajuste.

### Intervención causal limitada: pesos y estadísticas congelados

Para medir sensibilidad se sustituyó exclusivamente la característica **original** 15 por cero. En coordenadas normalizadas esto es `-muX(15)/sigmaX(15)`, no necesariamente cero. Se conservaron pesos, estadísticas y todas las otras entradas. Esta intervención **no es un modelo reentrenado válido para reemplazar al histórico**: sirve para aislar el efecto directo del canal.

La auditoría de los 49 bundles muestra:

- P20 a P480: ninguna conexión activa desde el canal 15.
- P490 y P500: una conexión activa.
- Fuente densa: las 12 conexiones están presentes; máximo peso absoluto aproximadamente `0.0479`.
- P500: sobrevive la conexión a la neurona 12, peso `-1.8656e-4`. Contribución RMS `1.8758e-4`, frente a `0.87349` de las otras entradas: razón `2.1474e-4`.

| Modelo | Rotación | E baseline | E tras intervención congelada |
|---|---:|---:|---:|
| Densa | 30° | `1.24476e-3` | `3.57207e-8` |
| Densa | 60° | `1.24736e-3` | `3.56590e-8` |
| Densa | 90° | `8.90387e-4` | `6.15366e-17` |
| Densa | 135° | `1.25056e-3` | `3.56439e-8` |
| Densa | 180° | `1.13060e-3` | `1.28465e-16` |
| Densa | 270° | `1.21062e-3` | `1.88871e-16` |
| P310 | 30° | `3.57395e-8` | `3.57395e-8` |
| P500 | 30° | `4.60205e-5` | `3.57125e-8` |
| P500 | 90° | `3.08335e-5` | `6.16116e-17` |

El cambio relativo de salida a 0° provocado por la intervención es `8.82680e-4` en la densa, **exactamente cero** en P310, y `3.24940e-5` en P500. La pérdida de equivariancia de la densa es medible y muy superior al suelo de precisión; en P500 es mucho menor. El residuo del orden de `3.6e-8` tras la intervención es compatible con la inferencia neuronal en precisión simple. Los ángulos de cuadrante no deben usarse por sí solos para estimar ese suelo.

![Error de equivariancia a 30 grados, con y sin el canal 15](../results/pnnn_phase_diagnostic/tp8c67afca_f6a5_423a_b86c_f1e6bfcf8106/phase_equivariance_channel15.png)

La figura usa la misma escala y no conecta modelos de distinta complejidad como si fueran una evolución temporal. Los símbolos se separan horizontalmente para mostrar la coincidencia en P310. También se guardan FIG y TikZ; el FIG se ha abierto correctamente.

### NMSE y límites de la conclusión causal

NMSE en dB a 0°, calculada con las mismas señales y filas:

| Modelo / intervención | Identificación | Complemento | Señal completa |
|---|---:|---:|---:|
| GMP-DOMP P310 histórico | -38.610857 | -38.403817 | -38.424040 |
| PN-IQ-GMP P310 histórico | -39.054649 | -38.864418 | -38.883032 |
| PNNN P310 histórico | -35.502860 | -35.482326 | -35.484371 |
| PNNN P310, canal anulado con pesos congelados | -35.502860 | -35.482326 | -35.484371 |
| Densa histórica | -35.451178 | -35.457236 | -35.456631 |
| Densa, canal anulado con pesos congelados | -35.446561 | -35.462126 | -35.460570 |
| P500 histórico | -35.756817 | -35.694609 | -35.700777 |
| P500, canal anulado con pesos congelados | -35.756810 | -35.694629 | -35.700794 |

En la densa, anular el canal **empeora** identificación aproximadamente 0.00462 dB y mejora el complemento 0.00489 dB; no es una mejora uniforme. En P500 las diferencias son del orden de `1e-5` dB. P310 no cambia en absoluto, no solo al redondear la tabla. La rotación de 30° altera la NMSE completa de la densa de -35.456631 a -35.443351 dB, mientras que la intervención la estabiliza al nivel de precisión numérica.

Se demuestra sensibilidad directa en densa/P500, pero **no se demuestra que explique la brecha de varios dB de P310**. Un posible efecto indirecto sobre la trayectoria del entrenamiento denso y la poda posterior queda pendiente del experimento desde cero. La mejora de equivariancia no demuestra por sí sola mejora de NMSE, y ninguna de estas medidas demuestra linealización RF.

## C. Alternativas y cambio experimental

### A: cero estructural exacto — implementada como opción explícita

`cfg.pnnn.zeroCurrentImaginaryTap = true` hace que `buildPhaseNormDataset`, solo en modo `full`, asigne `features(M+2,:) = 0`. No cambia las otras 83 posiciones, el target rotado, la arquitectura ni la dimensionalidad. La opción no existe en la configuración histórica por defecto: llamadas de cinco argumentos, modelos antiguos y valor explícito `false` conservan el cálculo anterior.

Se pasa esta misma opción desde entrenamiento, reconstrucción de características al reanudar y ONLINE. El modelo nuevo la almacena dentro de `featureConfig`. Con el entrenamiento nuevo, la estadística del canal es `muX=0`, `sigmaX=1` por la protección histórica de cero exacto; ONLINE usa precisamente esas estadísticas.

**No activar esta opción sobre un bundle histórico para presentarlo como modelo corregido.** El ejecutor experimental entrena su propia fuente densa y obtiene sus propias estadísticas. La intervención anterior solo diagnostica sensibilidad con pesos congelados.

Las posiciones de parámetros siguen siendo 1046 antes de podar y exactamente P según las máscaras después de podar. Sin embargo, un peso activo conectado a una entrada estructuralmente cero no influye en la función y no representa un grado de libertad efectivo. Puede conservar un valor inicial y sobrevivir a la poda por magnitud. No se fuerza su máscara a cero porque eso cambiaría el protocolo. En una futura comparación se deben comunicar tanto el presupuesto contado como las conexiones activas de esta columna; **P no equivale automáticamente a P grados de libertad efectivos**. El contador de FLOPs histórico también se conserva, sin introducir optimizaciones especiales del canal cero.

### B: tratamiento de sigma casi nula — analizada, no implementada

Una regla numéricamente defendible tendría que comparar sigma con una cota de error de la operación que produce cada característica, respetando sus unidades. Para este producto complejo la escala pertinente es `eps(double)*RMS(x)` y el residuo observado es compatible con ella. Para `abs(x)^p` la escala y propagación del error son diferentes. El cociente `std(feature)/RMS(feature)` aislado tampoco detecta un canal hecho enteramente de ruido de redondeo: aquí es aproximadamente uno.

Un suelo genérico de sigma sustituiría una identidad matemática exacta por una tolerancia elegida, podría afectar características legítimas de orden alto y requeriría almacenar/aplicar la misma regla ONLINE. Para este canal conocido, A es más sencilla y más justificable. No se ha añadido ningún umbral global ni cambiado el normalizador histórico.

## D. Experimento preparado y validación

### Comparación real P310: preparada, NO ejecutada

`runPNNNPhaseExperiment` reutiliza captura centrada, split y configuración del experimento guardado. Entrena una nueva PNNN desde cero con la opción activa, ejecutando selección interna de épocas, refit denso y poda/fine-tuning P310. Mantiene arquitectura, semillas, reglas de poda y presupuesto. **No congela las épocas a 1228/164**: mantiene el mismo criterio de selección, que puede escoger otras épocas con las nuevas características.

No repite DOMP ni modifica los modelos lineales: copia los dos modelos de control sin cambios en un nuevo bundle. Crea un directorio único bajo `results/pnnn_phase_experiment/` y guarda configuración/procedencia, fuente densa, punto P310, bundle ONLINE y `comparison.csv`. Las épocas y el tiempo se guardan en `bundle.phaseExperiment`; las métricas incluyen identificación, complemento, señal completa, parámetros y FLOPs. El diagnóstico posterior recoge pesos/máscaras, sensibilidad y rotaciones de la variante.

Coste orientativo: en el barrido guardado transcurrieron 8 min 16 s entre el checkpoint lineal fijo y la fuente densa, y 61 s entre los puntos P300 y P310. La prueba de una nueva fuente más un P310 tiene un orden de coste de **9–15 minutos** en condiciones semejantes; no es una garantía ni un benchmark dedicado, y las épocas seleccionadas pueden cambiar. Este coste se informa antes de lanzar el entrenamiento real.

Desde la raíz del repositorio, cuando se decida ejecutar la prueba:

```matlab
[modelFile, comparison] = runPNNNPhaseExperiment( ...
    'results/dpd_sweep/tp98e6f7ac_f3cd_44c7_9ce3_fd8619319748', 310);

report = diagnosePNNNPhaseNormalization(modelFile, ...
    'measurements/experiment20260429T134032_xy.mat');
```

No son comandos de barrido: entrenan una única fuente PNNN y un único presupuesto. El diagnóstico adicional solo infiere. No se ha reutilizado el complemento para ajustar variantes o hiperparámetros; la elección de cero estructural deriva de la identidad matemática y de la intervención causal, no de maximizar su NMSE.

### Pruebas realizadas

| Comprobación | Resultado |
|---|---|
| `run_pnnn_structural_zero_test` | PASS: captura sintética de 512 muestras, entrenamientos de 2 épocas; las otras 83 características y los targets no cambian |
| Opción ausente frente a `false` | Igualdad exacta de características; igualdad exacta de predicciones reales P310 |
| Variante sintética entrenada y guardada | `muX(15)=0`, `sigmaX(15)=1`, opción persistida; controles lineales idénticos, comparados incluyendo NaN en descriptores |
| Presupuesto y máscaras sintéticos | P20 y P310 exactos; P310 tiene 14 biases y 780 FLOPs; pesos fuera de máscara a cero |
| Entrenamiento experimental perfilado | Dos llamadas a `trainnet` para selección/refit; ninguna a `selectDOMPSupport` |
| Predicción directa frente a bundle recargado | Igualdad exacta |
| `run_dpd_online_fresh_session_test` en proceso MATLAB nuevo | PASS: predicciones exactas, sin DOMP, entrenamiento, poda ni fine-tuning |
| Rotaciones sintéticas 30/60/90/135/180/270° | `E < 1e-6`, tolerancia conservadora para la inferencia de red en precisión simple |
| P310 histórico, los tres modelos | Reproducción exacta (`isequal`) frente a las señales guardadas antes del cambio productivo |
| P340 histórico, los tres modelos | Reproducción exacta frente a `dpd_signals.mat` histórico |
| `run_dpd_model_reuse_test` | PASS |
| `run_pnnn_shared_dense_sweep_test` | PASS |
| `run_scientific_pipeline_smoke_test` | PASS |
| `run_dpd_parameter_sweep_test` adicional | Incompleto: detenido al alcanzar el límite de 180 s, durante el ajuste lineal de su segunda captura sintética. Las aserciones anteriores, incluidas reanudación completa/parcial y equivalencia de predicciones, habían terminado sin error. No se presenta como PASS global |
| `checkcode` de los nueve archivos MATLAB afectados/nuevos | Sin mensajes tras corregir una sugerencia de indexación en la auditoría |
| `git diff --check` | Sin errores |

La tolerancia `1e-6` solo se usa para equivariancia de la variante sintética; las comprobaciones de compatibilidad histórica y guardado/carga exigen **igualdad exacta**, no una tolerancia que pudiera ocultar cambios. Los tests sintéticos no predicen el rendimiento del entrenamiento real corregido.

El P340 comprobado es `results/dpd_offline/P340_20261011_210433_267_ed6696d7-1826-412c-b689-f7f61a53abfa/models.mat`, comparado con su inferencia histórica `online/inference_20261011_212149_467_46617fa4-db5e-477a-8571-1ca7278bcc17/dpd_signals.mat`.

### Artefactos de esta investigación

- [Diagnóstico P310, CSV/MAT y figura](../results/pnnn_phase_diagnostic/tp8c67afca_f6a5_423a_b86c_f1e6bfcf8106/): estadísticas de las 84 características, conexiones, rotaciones; `source_audit.mat`, `dense_rotations.csv` y `mask_audit.csv` con los 49 presupuestos.
- [Diagnóstico P500, CSV/MAT](../results/pnnn_phase_diagnostic/tp5011dea4_58a1_4f3f_8f45_11ca7f150046/): rotaciones 0/30/90°, conexiones y sensibilidad.
- [Experimento sintético validado, no modelo real P310](../results/pnnn_phase_experiment/tp177f330e_d021_4ef1_9ffa_b5ffab3797e7/): prueba de integración P20. El P310 sintético se probó en memoria partiendo de la misma fuente diminuta.

Estos directorios son nuevos. Los bundles/checkpoints/barridos históricos se han abierto solo para lectura. No se han cambiado sus parámetros, NMSE, soportes, lambdas, normalizaciones ni señales. Tampoco se han regenerado sus figuras. El WIP previo del repositorio se conserva; no hay commit ni push.

Archivos de código modificados en esta investigación: `toolbox/pnnn/buildPhaseNormDataset.m`, `toolbox/pnnn/prepare_pnnn_dense_source.m`, `toolbox/sweep/applyDPDModels.m`, `runDPDParameterSweep.m` y `run_parameter_sweep.m`. Salvo la opción de seis líneas en el constructor, los cambios son pasar `cfg.pnnn`/`featureConfig` a ese constructor. Archivos nuevos: `runPNNNPhaseExperiment.m`, `toolbox/pnnn/diagnosePNNNPhaseNormalization.m`, `toolbox/pnnn/auditPNNNPhaseSource.m`, `tests/run_pnnn_structural_zero_test.m` y este informe. Los otros cambios visibles en Git ya eran WIP previo y no forman parte de esta investigación.

## E. Conclusión científica

1. **Sí existe un problema numérico real:** el z-score transforma redondeo de una identidad cero en una entrada espuria de orden unidad.
2. **P310 no presenta efecto directo de ese canal en inferencia:** conexiones podadas, intervención exactamente nula, equivariancia en el suelo numérico. No explica de forma demostrada su peor NMSE.
3. **La densa sí pierde equivariancia de forma medible** (hasta `1.25e-3`); P500 muestra un efecto pequeño (`4.60e-5` a 30°). La intervención aislada reduce ambos al entorno de `3.6e-8`.
4. **Hay cambios de NMSE pequeños y no uniformes bajo intervención**, no una mejora demostrada del modelo reentrenado. El posible efecto sobre la trayectoria de entrenamiento/poda queda abierto.
5. **Conviene ejecutar primero la comparación única P310 preparada**, comunicando también conexiones estructuralmente inefectivas y épocas elegidas. No hay evidencia suficiente para justificar todavía repetir los 49 presupuestos. No se afirma ninguna mejora RF.
