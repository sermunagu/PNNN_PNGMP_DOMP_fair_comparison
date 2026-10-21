# Hoja de ruta — DPD con GMP-DOMP, PN-IQ-GMP y PNNN

Actualización contrastada el **21/10/2026**. Objetivo orientativo: Transactions
en enero de 2027; no es un plazo editorial confirmado. Repositorio local:
`C:\Sergi\Investigacion\Códigos\PNNN_PNGMP_DOMP_fair_comparison`.
El WIP local es la fuente de verdad, no GitHub. Este plan no autoriza nuevos
entrenamientos extensos, transmisiones RF, commits ni pushes.

La versión anterior completa, incluidos su registro histórico y TODO, queda
preservada en [plan_before_closure.md](../results/gmp_phase_study/20261020_205246/plan_before_closure.md).
Se corrigen hitos obsoletos: el sweep real, P310, la ablación en tres P y el
experimento PNNN de cero estructural **ya terminaron**. No hay que repetirlos.

## 0. Estados, evidencia y prioridades

- **COMPLETADO**: evidencia suficiente para el alcance indicado, no para uno mayor.
- **PARCIAL**: existe evidencia, pero falta una condición importante.
- **PENDIENTE OFFLINE**: trabajo posible sin nuevas transmisiones.
- **PENDIENTE RF**: requiere nuevas medidas o confirmación del banco con el tutor.
- **OPCIONAL**: no es requisito para defender la evidencia actual.

Prioridad 1: argumento y redacción PN-IQ con evidencia existente.
Prioridad 2: confirmación independiente/RF. Prioridad 3: extensiones.
Un comentario respondido en un informe no se considera incorporado al paper:
evidencia y redacción se clasifican por separado.

| Clave | Evidencia local |
|---|---|
| S | `results/dpd_sweep/tp98e6f7ac_f3cd_44c7_9ce3_fd8619319748/`: 49 presupuestos, 147 filas, bundles, configuración, split y checkpoints |
| A | [Ablación PN-IQ](PNIQ_GMP_ABLATION_2026-10-12.md), cierre J; `results/pniq_ablation/summary_20261013_164459/` y MAT P100/P310/P500 |
| N | [Diagnóstico PNNN](PNNN_PHASE_NORMALIZATION_DIAGNOSIS_2026-10-12.md); resultado real posterior en `results/pnnn_phase_experiment/tp1bfa4378_d472_497b_87d7_b5c638c34392/comparison.csv` y [guía actual](GUIA_REPOSITORIO_DPD.md) |
| E | [Eficiencia NMSE–FLOPs](DPD_EFFICIENCY_ANALYSIS.md); `results/dpd_efficiency/20261020_204728/` |
| Q | [Equivariancia y restricción GMP](GMP_PHASE_EQUIVARIANCE.md); `results/gmp_phase_study/20261020_205246/` |
| R | Consulta R2026a real: `regularization.mat` y `regularization_matlab.log` dentro de Q |
| O | [Flujo DPD](DPD_ILC.md), [guía](GUIA_REPOSITORIO_DPD.md), código OFFLINE/ONLINE y verificaciones históricas documentadas |

N contiene secciones anteriores al experimento real. El estado actualizado
se toma del CSV y la guía, no de la antigua frase «preparada, no ejecutada».

## 1. Alcance y contrato científico

En `measurements/experiment20260429T134032_xy.mat`, el mapeo
`xy_forward` aprende **u → x_ILC**: x es referencia deseada, y predistorsión
objetivo. La salida inferida de un DPD no es la salida medida del PA.
El archivo `_forward_xy.mat` del modelado PA tiene otra función física.
Confirmar siempre la señal realmente enviada al equipo, incluido posible CFR.

Se conservan 491520 muestras, 49151 de identificación y 442369 de evaluación.
Identificación/evaluación/full se reportan por separado. La señal completa
incluye entrenamiento. El complemento **no es una captura independiente**:
comparte adquisición, DC global, selección del segmento por el máximo de y y
memoria periódica; además se usa para seleccionar P. Ninguna cifra offline
demuestra por sí sola ACLR/ACPR/EVM o linealización RF.

Cada captura nueva que se vaya a identificar necesita entrenamiento nuevo
en carpeta propia, sin reutilizar automáticamente soportes, pesos ni
normalizaciones de otra captura. Sí se reutilizan etapas/checkpoints dentro
de una misma ejecución. ONLINE significa inferencia sin entrenamiento, no
garantía de procesamiento RF en tiempo real.

Se preservan PA P340, DPD recuperado P340, sweep DPD P310, los 15 modelos de
ablación y variantes científicas. Los restringidos no sustituyen el baseline.

## 2. Fase A — OFFLINE/ONLINE

| Tarea original | Estado | Evidencia | Qué falta / siguiente acción | Prioridad |
|---|---|---|---|---:|
| A1: archivo absoluto/relativo x/y, conservar fs | COMPLETADO | `DPD/trainDPD.m`, O | Semántica por campaña, no deducirla del nombre | 2 |
| A2: tres modelos nuevos; densa, poda y fine-tuning PNNN a P | COMPLETADO | `trainDPD`, S, pruebas sintéticas previas | No reentrenar por este cierre | 1 |
| A3: bundle autocontenido, métricas/configuración/procedencia y no sobrescritura | COMPLETADO | Bundles S y O | Mantener carpeta propia por ejecución | 1 |
| A3: ONLINE por archivo x sin y, mediante `applyDPDModels` | COMPLETADO | `DPD/runDPDOnline.m`, O | Ninguno para entrada compatible | 1 |
| A4: dos capturas sintéticas independientes, P exacto, tres salidas complejas, recarga y cero entrenamiento ONLINE | COMPLETADO | Verificaciones previas O, no repetidas ahora | No confundir test sintético con RF | 1 |
| A4: P340 y seis predictores históricos; sesión nueva | COMPLETADO | O, igualdad exacta reportada; archivos preservados | Ninguno en este cierre | 1 |
| Primer entrenamiento/barrido real autorizado | COMPLETADO | S: 49 puntos por arquitectura | Nueva captura sería otra campaña, no un TODO de software | 1 |
| Contrato de medidas: x/y, fs, DC, periodicidad, amplitud y sincronización | PENDIENTE RF | Convención actual documentada, no contrato universal | Confirmar con tutor por campaña | 2 |
| Límites DAC, CFR/PAPR y señal realmente enviada | PENDIENTE RF | No validación de banco equivalente en E/Q | Revisar antes de transmitir | 2 |
| NMSE opcional de otra entrada ONLINE con objetivo disponible | OPCIONAL | Inferencia explícita disponible | Evaluación adicional, no métrica RF | 3 |

Dos comandos del flujo normal, desde la raíz y con DPD en el path:

```matlab
modelFile = trainDPD('medidas_nuevas.mat', 310);
runDPDOnline(modelFile, 'entrada_nueva.mat');
```

El primero **entrena realmente** y no se ha ejecutado en este cierre.
El segundo usa x y guarda `dpd(k).yvalmod` / `dpd(k).modeltype` en archivo
nuevo. Confirmar el contrato físico antes de RF.

## 3. Fase B — Validación científica

| Tarea original | Estado | Evidencia | Qué falta | Prioridad |
|---|---|---|---|---:|
| B1: protocolo histórico, identificación/complemento/full y mismo NMSE/índices para las tres familias | COMPLETADO | S, A, E, Q | Ninguno para reportar esta captura | 1 |
| B1: validación interna PNNN separada del complemento | COMPLETADO | `trainDensePNNN`, S | Explicar selección provisional/refit en paper | 1 |
| B1: generalización independiente | PARCIAL | Complemento disponible, no test externo | Captura independiente (RF), o protocolo predefinido con preprocesado/guardas y nueva validación autorizada | 2 |
| B2: referencia sin DPD, ILC y tres DPD en mismo montaje/estado | PENDIENTE RF | S aporta modelos, no salidas RF comparables | Adquirir referencias y salidas PA reales | 2 |
| B2: entrada enviada, salida, ganancias, potencia, fs, CFR/DAC, alineación | PENDIENTE RF | `rf_measurements.csv` sin asociaciones RF | Completar por medida sin sustituir inferencias | 2 |
| B2: NMSE RF, ACPR/ACLR, espectros medidos y EVM aplicable | PENDIENTE RF | Sin resultados RF de esta campaña en S/E/Q | Definir bandas/referencia/procedimiento y medir | 2 |
| B2: comparabilidad waveform/potencia/configuración/adquisición | PENDIENTE RF | Requisitos descritos, no campaña terminada | Protocolo acordado y repetición controlada | 2 |
| B3: canal PNNN 15 casi nulo y normalización | COMPLETADO | N real P310: -35.482326 → -35.525255 dB; mismas 310/780 | Ganancia 0.042929 dB; no justifica otro sweep neuronal | 1 |
| B3: normalización, condición y sensibilidad | COMPLETADO | A: SVD de 15 modelos y perturbación dirigida; Q: seis restringidos | No prueba robustez universal ni en cuantización | 1 |
| B3: características, retardos, operaciones especiales y poda ideal | COMPLETADO | E y desgloses guardados | Hardware solo opcional; conservar cifras históricas | 1 |
| B3: interpretación física de memoria/amplitud/fase | PARCIAL | E por amplitud; Q algebraico/experimental | Causalidad física requiere medidas controladas | 2 |

La variante PNNN real seleccionó 1221/163 épocas frente a 1228/164 del
baseline. No se entrenó PNNN en E/Q ni en este cierre. La variante estructural
no reemplaza la normalización histórica.

## 4. Fase C — Barrido DPD: completado offline

| Tarea original | Estado | Evidencia | Qué falta | Prioridad |
|---|---|---|---|---:|
| Ejecutor, una ruta DOMP por arquitectura/captura y prefijos | COMPLETADO | `DPD/runDPDSweep.m`, S | Ninguno | 1 |
| Densa definitiva común; poda/fine-tuning independientes por P | COMPLETADO | S y verificaciones O | Ningún entrenamiento adicional solicitado | 1 |
| Checkpoints y reanudación explícita sin repetir etapas | COMPLETADO | S, O | Mantener compatibilidad | 1 |
| Bundles compactos/remapeo PN-IQ, P, FLOPs, tres NMSE, rangos, sparsity y épocas | COMPLETADO | 49 bundles S y tablas | Ninguno para campaña offline | 1 |
| Figuras complejidad/sparsity/calidad | COMPLETADO | S, E | Selección editorial | 1 |
| Ejecutar campaña real 20:10:500 | COMPLETADO | 147 filas S | No repetir | 1 |
| Seleccionar punto con malla completa | COMPLETADO | DPD P310 y PA legacy P340; selectores comprobados | P310 es punto de desarrollo, no óptimo RF demostrado | 1 |
| Puntos RF y calidad RF vs P/FLOPs | PENDIENTE RF | P100/P310/P500 son candidatos documentados | Acordar y medir, no transmitir 49 por defecto | 2 |

Se conserva el selector corregido: la malla piloto dispersa no autoriza
ventanas ficticiamente completas. No se altera P310.

## 5. Fase D — Eficiencia NMSE–FLOPs

| Tarea original | Estado | Evidencia | Qué falta | Prioridad |
|---|---|---|---|---:|
| Convención FLOPs y costes excluidos | COMPLETADO | E y contadores actuales | Exponerla en paper, no llamarla tiempo real | 1 |
| Sparse ideal vs ejecución densa con máscaras | COMPLETADO | PNNN P310: 780 ideal, 2252 núcleo denso bajo contador | Sin benchmark físico | 1 |
| Igual P y tablas sin reentrenar | COMPLETADO | S/E | Redacción | 1 |
| Pareto por familia/global y dominados | COMPLETADO | E: global 49 PN-IQ + GMP P20/P30/P40; 95 dominados | No extrapolar a otras capturas | 1 |
| Igual techo de FLOPs con modelos discretos | COMPLETADO | E: 145 techos, 435 filas; PN-IQ gana desde F=245 en techos observados | Ningún modelo interpolado ni RF inferido | 1 |
| Igual NMSE objetivo | COMPLETADO | E: -38 dB exige GMP 812, PN-IQ 486; PNNN no alcanza | Umbrales RF aparte | 1 |
| Error por amplitud, energía y muestras por bin | COMPLETADO | E: 20 cuantiles P310, 22118–22119 muestras/bin | No inferir causa física de la correlación | 1 |
| Igual calidad RF y distorsión fuera de banda | PENDIENTE RF | Sin asociación RF S | ACPR/ACLR/EVM medidos | 2 |
| Latencia/energía/recursos de implementación | OPCIONAL | FLOPs no los prueban | Benchmark si se formula esa afirmación | 3 |

PN-IQ P310 gana 0.460601 dB a GMP con 38.63% menos FLOPs. A techo 780:
PN-IQ P220 obtiene -38.639739 dB, GMP P140 -37.993928 y PNNN P310
-35.482326. No se exige gastar exactamente el techo. En las siete primeras
bandas de amplitud (35% de muestras) GMP supera a PN-IQ; PN-IQ gana en las
otras trece. El error se concentra en amplitudes altas, sin demostrar memoria
térmica, load modulation u otra causa física.

No se suman indiscriminadamente raíces/divisiones/exponenciales a FLOPs.
El coste ideal PNNN omite pesos cero y no incluye todo el z-score/activación;
los campos genéricos del contador conservan el núcleo denso. E detalla estos
límites. No se corrigen retrospectivamente 780, 1023 ni 1667.

## 6. Fase E — Ablación y equivariancia

| Tarea original | Estado | Evidencia | Qué falta | Prioridad |
|---|---|---|---|---:|
| Definir GMP/PN-GMP/IQ/PN-IQ sin variantes artificiales | COMPLETADO | A: rotación unitaria compleja equivale al LS original; IQ nativo/común definidos | No presentar PN-GMP equivalente como factor independiente | 1 |
| Control P, diccionario y soporte común | COMPLETADO | A: 15 modelos P100/P310/P500, mapas I/Q válidos | Redactar sesgo del soporte seleccionado para IQ | 1 |
| Separar PN de «más grados de libertad» | PARCIAL | Mismo soporte/P: ganancias 0.495203/0.448105/0.596964 dB | No es descomposición causal total ni evidencia universal | 1 |
| Interacción factorial aditiva PN × IQ | OPCIONAL | Equivalencia compleja impide interpretar cuatro nombres como factores independientes | Definir primero estimando/controles; no afirmar interacción medida | 3 |
| Adaptar probe widely-linear PA a DPD | OPCIONAL | Probe PA conservado; A aporta controles DPD | Solo para una pregunta distinta bien definida | 3 |
| Comparar ablación a igual P y costes | COMPLETADO | A tablas P/FLOPs; E familias principales | Tres P de controles no equivalen a curva densa de cada control | 1 |
| Clasificar auxiliares y diagnosticar congelados | COMPLETADO | Q: 673 fuentes, 671 grado fase 1, conj(x) grado -1, abs(x) grado 0 | Ninguno para esta población | 1 |
| Nuevo DOMP restringido y LS P100/P310/P500 | COMPLETADO | Q: dos rutas nuevas, seis modelos, soportes/SVD separados | No sustituir históricos | 1 |
| Compromiso simetría/precisión | COMPLETADO | Q: defecto fase ~0.003 → ~1e-11, cinco de seis NMSE empeoran | Externa/RF si se propone usar la restricción | 2 |

PN-IQ restringido gana a GMP en P310/P500 (0.384481/0.539807 dB) y pierde en
P100 (0.265394 dB). Los auxiliares no explican toda la ventaja PN-IQ.
La identidad de equivariancia no es novedosa por sí misma ni está garantizada
para PN-IQ en muestras exactamente nulas con taps no nulos. En esta captura
no hay ceros exactos; Q documenta excepción y redondeo. No se demuestra que
los auxiliares históricos sean físicamente erróneos.

## 7. Fase F — Reviewers: evidencia frente a respuesta en paper

Solo se atribuyen a reviewers las **cinco observaciones recopiladas en la
versión anterior de este plan**. No se ha verificado aquí una nueva carta
original ni una respuesta enviada. POMP/Volterra corresponde a otros
manuscritos: no se incorporan comentarios supuestos.

| Origen / observación recopilada | Estado | Evidencia disponible | Qué falta | Prioridad |
|---|---|---|---|---:|
| R1: regularización PNNN, falta valor en tabla | PARCIAL | **Consulta técnica COMPLETADA** R2026a: L2 densa=1e-4; peso=1, sesgo=0; fine-tuning sin L2 añadida | Incorporar valor/distinción de etapas a manuscrito, tabla y respuesta | 1 |
| R1: añadir ACPR | PENDIENTE RF | Sin medidas comparables asociadas | Medir/publicar ACPR/ACLR real con bandas definidas | 2 |
| R1: distorsión fuera de banda vs complejidad PNNN | PARCIAL | Sweep/coste offline S/E completados | RF/espectros medidos, no inferir distorsión PA del error ILC | 2 |
| R2: identificar ~10% y evaluar sobre 100% | PARCIAL | Tres métricas y complemento comunes reportados | Límites, captura externa y protocolo de confirmación; trasladar al paper | 1/2 |
| R2: ganancia PN frente a grados de libertad | PARCIAL | A a igual P/soporte, equivalencia unitaria y Q | Respuesta sin atribución causal total ni factorial ficticio | 1 |

### Regularización efectiva, comprobada directamente

`trainingOptions("adam")` en **MATLAB R2026a** devuelve
`L2Regularization=1e-4`. `fullyConnectedLayer(2)` devuelve
`WeightL2Factor=1`, `BiasL2Factor=0`. Ambas llamadas a `trainingOptions`
en `toolbox/pnnn/trainDensePNNN.m` omiten una sobrescritura L2 y fc1/fcOut
usan los factores por defecto. `trainnet` recibe `"mse"`. Por tanto, el
entrenamiento denso provisional y definitivo utiliza MSE con regularización
L2 sobre pesos, no sobre sesgos: **no es L2=0**. La red conserva fc1,
sigmoide y fcOut. No se modificó arquitectura ni métrica.

`toolbox/pnnn/pruneAndFineTunePNNN.m` calcula
`mean((dlY-dlT).^2,"all")`, deriva ese loss, aplica máscaras y usa
`adamupdate`; no añade L1/L2 ni weight decay. Los factores L2 de las capas
no añaden por sí solos penalización al loss del bucle personalizado.
Poda por magnitud y selección de épocas son mecanismos distintos de L1/L2.
La validación interna selecciona duración; el refit usa épocas elegidas,
no validación sobre el complemento. No se entrenó para esta consulta.

Alcance: código actual + defaults efectivos R2026a. No es reconstrucción
forense de opciones no serializadas en cada entrenamiento histórico.
R conserva valores y consulta; no se cambió checkpoint ni hiperparámetro.

Condicionamiento, normalización, equidad y FLOPs son requisitos técnicos
revisados en B3/D/E, sin inventar atribuciones adicionales a R1/R2.

## 8. Fase G — Preparación del artículo

| Tarea original | Estado | Evidencia | Qué falta | Prioridad |
|---|---|---|---|---:|
| Pregunta central/contribución PN-IQ | PARCIAL | Propuesta abajo, E/A/Q | Integrar en manuscrito sin novedad matemática no demostrada | 1 |
| Separar PA/imitación ILC/linealización PA | PARCIAL | O y este plan | Revisar cada tabla/figura del paper | 1 |
| Dataset, ILC, variables y preprocesado | PARCIAL | Captura/split/configuración trazables | Texto métodos; calibración/sincronización banco con tutor | 1/2 |
| Comparación justa P/FLOPs/calidad | PENDIENTE OFFLINE | E terminado | Integrar resultados y límites del contador en paper | 1 |
| Ablación/grados de libertad | PENDIENTE OFFLINE | A/Q terminados | Incorporar controles, negativos y límites | 1 |
| Generalización/evaluación no entrenada | PARCIAL | Complemento; sin confirmación externa | Nueva captura/protocolo predefinido | 2 |
| ACPR/ACLR, EVM y espectros medidos | PENDIENTE RF | No cerrado | Adquisición/análisis RF | 2 |
| Incertidumbres memoria/potencia/normalización/hardware | PARCIAL | Límites E/A/Q documentados | Variación entre capturas/potencias; evitar causalidad no medida | 1/2 |
| Trazabilidad tablas/figuras/scripts/modelos | COMPLETADO | S/A/E/Q, scripts copiados y MAT/CSV | Conservar correspondencia al seleccionar material | 1 |

### Propuesta científica preliminar (aún no incorporada al paper)

Estudiamos PN-IQ-GMP como aproximador DPD eficiente de una predistorsión ILC,
comparándolo con GMP-DOMP complejo y PNNN sparse bajo protocolo común.
En esta captura, PN-IQ mejora el compromiso precisión/coste a presupuestos
moderados y altos: P310 aporta 0.4606 dB frente a GMP con 38.6% menos FLOPs
analíticos; a igual techo discreto también domina en la región observada
F≥245. No gana en todo P ni en toda región de amplitud.

La ablación a igual P y con los mismos identificadores de características
respalda un efecto de representación no reducible a rotar unitariamente el
LS complejo. Candidatos y selección siguen importando: no atribuimos toda
la ganancia a PN aislada. La restricción de auxiliares cuantifica el compromiso
equivariancia/precisión sin eliminar la ventaja PN-IQ en P310/P500.
Es evidencia de diseño, no innovación algebraica por el mero uso de PN.

La contribución más sólida hoy es **experimental offline y comparativa**.
Sería débil presentar la identidad de equivariancia como novedad o afirmar
linealización RF solo por imitar x_ILC. RF deberá confirmar ACPR/ACLR/EVM/NMSE
del PA con waveform/potencia/montaje comparables; otra captura debe confirmar
generalización. Esa contribución RF todavía no está establecida.

## 9. Decisiones abiertas: actualización de D1–D10

| ID | Estado | Decisión / pendiente |
|---|---|---|
| D1, x/y en nuevas medidas | PARCIAL | Contrato software disponible; significado físico/fs por campaña |
| D2, entrada ONLINE independiente | COMPLETADO | Archivo explícito, x, sin y obligatoria |
| D3, reentrenar misma captura | COMPLETADO | Petición explícita, ejecución nueva sin sobrescribir; distinto de recuperación |
| D4, NMSE principal | COMPLETADO | Complemento para desarrollo/selección; identificación/full secundarios; no test externo |
| D5, x frente a xCFR en PA | PENDIENTE RF | Confirmar señal enviada y adquisición con tutor |
| D6, cambio PA/nueva identificación | COMPLETADO | Nueva captura a identificar: entrenamiento nuevo; reutilización explícita para generalización |
| D7, FLOPs igualados | COMPLETADO | Techos discretos, mejor modelo factible; no interpolación como resultado |
| D8, cuatro variantes PN/IQ | COMPLETADO | Controles matemáticos A; no factorial aditivo medido |
| D9, metadatos RF | PENDIENTE RF | Potencia/ganancias/fs/alineación/CFR/DAC/bandas/referencias |
| D10, canal PNNN casi nulo | COMPLETADO | Diagnóstico y variante real P310 guardados; baseline intacto |

## 10. Próximos pasos priorizados

| Tarea | Estado | Evidencia | Qué falta | Prioridad |
|---|---|---|---|---:|
| Métodos, tabla eficiencia y respuesta limitada R1/R2 | PENDIENTE OFFLINE | E/A/Q/R y S | Redacción editorial/integración manuscrito | 1 |
| Seleccionar figuras no redundantes y resultados negativos | PENDIENTE OFFLINE | Cinco figuras E/Q y ablación existente | Selección, no nuevos cálculos | 1 |
| Predefinir validación externa/puntos RF | PARCIAL | P310, extremos P100/P500 conocidos | Acuerdo tutor: referencia/captura/potencia/métricas/límites | 2 |
| Nueva captura y linealización comparativa | PENDIENTE RF | Modelos/exportaciones disponibles | Laboratorio autorizado, calibración y datos | 2 |
| Analizar RF con trazabilidad/variación experimental | PENDIENTE RF | Plantilla `rf_measurements.csv` | Asociar medidas y calcular métricas reales | 2 |
| Precisión reducida/hardware/perturbación de matriz | OPCIONAL | Condiciones ~1e8 aconsejan prudencia | Solo para afirmaciones específicas de implementación | 3 |
| Más arquitecturas/podas PNNN u otro sweep de 49 redes | OPCIONAL | Baseline y variante P310 existen | No recomendado para este objetivo ni autorizado aquí | 3 |

## 11. Cierre de ejecución y conservación

El 20/10 se terminó E y se guardaron todos los cálculos Q. El 21/10 se
inspeccionó MATLAB PID32796, CPU, log y checkpoints: la exportación llevaba
más de 13 h sin nuevos archivos. Se cerró solo ese proceso identificado;
no se lanzó otro hasta confirmar su salida. Se completaron el TikZ de NMSE y
la figura de equivariancia P310 desde MAT/FIG, sin repetir cálculos.
Ambos gráficos Q se inspeccionaron visualmente; `matlab2tikz` avisa de que
`legend('Location','best')` usa una posición por defecto en TikZ, sin afectar datos.

Pasaron: tablas E iguales a S, splits E/Q idénticos, muestras por bins,
seis presupuestos activos P exactos/rangos P/2, longitudes 491520,
selectores DPD P310/PA legacy P340 y apertura de los cinco FIG. El primer
chequeo se detuvo después en tres avisos `checkcode` por supresiones AGROW
obsoletas; se retiraron solo esos comentarios, sin repetir pruebas superadas.
El `checkcode` final pasó sin mensajes en ambos ejecutores y en los scripts
de recuperación/consulta; los logs se conservan en Q.
R2026a confirmó L2=1e-4. El inventario de 511 archivos históricos
(medidas/resultados fuera de E/Q) coincide en rutas, tamaños y fechas con
el previo: control de no escritura, no nueva comparación criptográfica
ni repetición de inferencias.

Se mantienen las verificaciones exactas históricas de seis predictores
P310/P340 y quince de ablación documentadas, sin reejecutarlas.
No se modificaron algoritmos compartidos, normalizaciones, soportes,
coeficientes, lambdas, métricas ni checkpoints históricos. No hubo nuevo
entrenamiento PNNN, sweep completo, RF, commit ni push.

## 12. Referencias operativas actuales

Desde la raíz: `addpath('PA','DPD','PA/experiments','DPD/experiments')`.

- OFFLINE/ONLINE: `DPD/trainDPD.m`, `DPD/runDPDOnline.m`, `DPD/exampleDPD.m`, `toolbox/sweep/applyDPDModels.m`.
- Sweep DPD: `DPD/runDPDSweep.m`, reanudación por carpeta explícita; no ejecutar para este cierre.
- PA: `PA/runPA.m`, `PA/runPASweep.m`, `PA/comparePA.m`. `comparePA(P)` sin sweep puede iniciarlo completo: no usar como consulta.
- Experimentos: `DPD/experiments/runPNIQAblation.m`, `runPNNNPhaseExperiment.m`, `analyzeDPDEfficiency.m`, `runGMPPhaseStudy.m`.
- Configuración: `config/getComparisonConfig.m`; lineales `fitLinearModels`, `fitComplexGMP`, `fitPNIQGMP`, `fitFixedRidge`.
- Nombres anteriores se mantienen en los registros históricos como procedencia, no como comandos operativos.
