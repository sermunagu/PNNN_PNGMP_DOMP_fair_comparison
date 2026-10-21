# Ablación matemática y experimental de PN-IQ-GMP

> Organización actual (14/10/2026): ejecutores PA en `PA/` y DPD en `DPD/`;
> `runPNIQAblation` y `runPNNNPhaseExperiment` están en `DPD/experiments/`.
> Antes de usar ejemplos: `addpath('PA','DPD','PA/experiments','DPD/experiments')`.
> Las narraciones históricas conservan sus nombres; el [mapeo actual](GUIA_REPOSITORIO_DPD.md#entradas-reorganizadas-14102026)
> indica sus sustitutos, sin cambiar medidas, modelos ni resultados.

## A. Objetivo y alcance

Separar, hasta donde permite una comparación de soporte limitado, la representación compleja, la representación real I/Q, la normalización de fase y la selección DOMP. No se propone un factorial 2×2: GMP y PN-GMP complejos son algebraicamente equivalentes bajo rotación común de filas.

Se utiliza exclusivamente DPD por imitación del ILC: `x=u`, `y=x_ILC`, captura `measurements/experiment20260429T134032_xy.mat`, barrido `results/dpd_sweep/tp98e6f7ac_f3cd_44c7_9ce3_fd8619319748/`. No se usa la captura `forward_xy.mat` del PA. P310 es el punto operativo seleccionado previamente; no se vuelve a seleccionar P con esta ablación. No se entrena ni estudia PNNN en este bloque.

## B. Formulaciones y equivalencias

### B1. Rotación común: un control algebraico, no otro modelo

Para `r=conj(x)/abs(x)` cuando `x~=0`, y `r=1` cuando `x=0`, sea `D=diag(r)`. En aritmética exacta `D^H D=I`, por lo que

\[
\|Dy-DUc\|_2^2=\|y-Uc\|_2^2,\qquad D^HDUc=Uc.
\]

El pico y la norma 2 de cada columna también se conservan. El LS complejo con el mismo soporte y el mismo criterio de rango tiene la misma predicción restaurada. Si hay deficiencia de rango, el conjunto de minimizadores y el criterio de mínima norma de coeficientes tampoco cambian bajo esta transformación unitaria de filas.

DOMP tiene la misma invariancia en aritmética exacta: si `Z_PN=DZ` y `e_PN=De`, se conservan las normas de columnas, `Z_PN^H e_PN=Z^H e`, las puntuaciones y las proyecciones de Gram–Schmidt. El reajuste LS y el residuo se transforman de la misma forma. Por inducción, se conserva el orden de selección salvo la resolución de empates. En coma flotante, columnas mal condicionadas, empates o umbrales de rango pueden producir diferencias; no se presume identidad universal bit a bit de rutas reales largas.

La prueba MATLAB pequeña usa 128 filas y 12 columnas complejas: error relativo de la predicción restaurada `6.71e-16`, y los mismos seis índices DOMP. Su tolerancia `1e-12` es conservadora para esta matriz bien condicionada en doble precisión. No se mantiene un pipeline PN-GMP duplicado ni se repite DOMP complejo sobre la captura real.

Esta invariancia **no exige** equivariancia física de cada columna ante `x -> exp(j*theta)*x`: es una transformación de filas de una matriz ya construida, no una reconstrucción del diccionario con otra entrada.

### B2. Complejo acoplado frente a I/Q libre

Escribiendo `U=R+jJ`, `c=u+jv`, GMP impone:

\[
\hat y_I=Ru-Jv,\qquad \hat y_Q=Rv+Ju.
\]

Las dos componentes de cada columna están acopladas por los mismos dos coeficientes reales. IQ libre usa `A=[R,J]` y `yhat=Aa+jAb`, sin esos acoplamientos entre salidas. El soporte real se selecciona conjuntamente para ambas salidas; **no son dos soportes DOMP independientes**. Para una columna real, la puntuación con objetivo complejo equivale a la norma conjunta de sus correlaciones con I y Q.

Si están presentes todas las componentes, separando `alpha=a+jb=[alpha_R;alpha_J]`,

\[
A\alpha=Uc+\overline U d,\quad
c=(\alpha_R-j\alpha_J)/2,\quad d=(\alpha_R+j\alpha_J)/2.
\]

Es una representación widely-linear sobre el diccionario. Pero seleccionar 155 componentes reales y ajustar dos coeficientes por componente **no equivale** a seleccionar 155 átomos complejos libres de `[U,conj(U)]`: una componente real aislada liga sus ramas directa y conjugada. La familia I/Q completa contiene al GMP complejo completo, pero los modelos dispersos con **el mismo P** no forman necesariamente familias anidadas.

El probe existente `run_widely_linear_gmp_probe` no sirve como resultado de esta ablación: utiliza el mapeo/captura PA por defecto, P340, selección interna de lambda Ridge y átomos complejos `[U,conj(U)]`. Se ha leído, pero no ejecutado ni trasladado como si fuera el mismo experimento.

### B3. PN-IQ: una clase funcional diferente

Con `V=DU`, PN-IQ ajusta `[Re(V),Im(V)]` al objetivo `Dy` y restaura `D^H`. Su forma ampliamente lineal equivalente es

\[
\hat y=Uc+\overline r^{\,2}\overline U d.
\]

La rama conjugada incorpora una modulación dependiente de la fase instantánea de la entrada. Por eso la invariancia LS de B1 **no** demuestra equivalencia entre IQ y PN-IQ libres: descomponer en partes reales y permitir coeficientes independientes cambia las restricciones del problema.

Para un monomio de grado de fase `q=numel(X)-numel(Xconj)`, bajo rotación global sus ramas son:

| Formulación | Grados de fase de las ramas |
|---|---|
| GMP complejo | `q` |
| IQ sin PN | `q` y `-q` |
| PN-IQ | `q` y `2-q` |

Así, los términos GMP canónicos (`q=1`) conservan la equivariancia en PN-IQ aun con coeficientes I/Q libres; en IQ sin PN la rama conjugada puede romperla. Los auxiliares `conj(x)` (`q=-1`) y `abs(x)` (`q=0`) no satisfacen esa condición. La identidad se entiende en muestras con tap actual no nulo; la convención en cero se mantiene sin cambiar el código histórico.

## C. Diseño experimental y contabilidad

### Población, muestras y ajuste

La población histórica contiene 143 términos principales, 264 de envolvente retardada y 264 adelantada, más `conj(x)` y `abs(x)`: 673 tras eliminar el `x` duplicado. `modelsConfig` fija orden 13, lags de portadora 0:10 y desplazamientos de envolvente 1:2; `regressorManager` conserva su orden y canonicalización.

PN-IQ elimina 61 columnas Q de términos canónicos con portadora de lag cero: 13 principales, 24 retardadas y 24 adelantadas. La factorización analítica de `buildPNIQFeatures` fuerza esos ceros; no se utiliza esa función para fabricar IQ sin PN.

| Espacio candidato | Número | Observación |
|---|---:|---|
| GMP complejo | 673 | Coeficientes complejos |
| IQ nativo | 1345 | `[Re(U),Im(U)]`, eliminando solo `Im(abs(x))=0` |
| PN-IQ nativo | 1285 | `2*673-61` |
| Identificadores admisibles comunes | 1284 | Intersección por regresor fuente y componente I/Q |

IQ nativo conserva candidatos exactamente redundantes, como `Re(conj(x))=Re(x)` y `Im(conj(x))=-Im(x)`; DOMP trata su dependencia al ortogonalizar. No se aplica a IQ la eliminación de 61 columnas propia de PN.

Se reutilizan las 491520 muestras centradas y el split guardados en `experiment.mat`: 49151 de identificación, 442369 de complemento. Los taps son periódicos. No se estima otra ganancia, escala, semilla o partición. La captura completa participa en el centrado DC y en la selección histórica del segmento de identificación; esto limita la independencia de la evaluación.

Los ajustes nuevos utilizan el DOMP existente y el LS principal vigente: **pico por columna y pico de salida**, `lsqminnorm` sin Ridge, tolerancia `max(size(A))*eps(norm(A,2))` sobre la matriz normalizada, y desnormalización de coeficientes antes de predecir. No se usa la antigua definición de norma 2 para el ajuste final. DOMP conserva su normalización interna a norma 2 y su reajuste QR/fallback LS; son etapas distintas.

### Conjunto mínimo de comparaciones

1. GMP-DOMP histórico, sin reajuste.
2. PN-IQ histórico, sin reajuste.
3. IQ nativo, DOMP propio y 155 características reales.
4. IQ restringido a candidatos comunes, DOMP propio y 155 características.
5. PN-IQ sobre exactamente los 155 identificadores seleccionados en 4, sin repetir DOMP.

El último control favorece potencialmente a IQ porque su soporte se eligió allí. No es un soporte neutral. La comparación 3–4 mide el cambio de candidatos **junto con su efecto en la ruta DOMP**, mientras que 4–5 mantiene el soporte fijo y cambia representación y objetivo. Comparar 2–5 mezcla selección y la disponibilidad del candidato PN exclusivo; no se atribuye su diferencia exclusivamente a DOMP.

El traslado inverso del soporte PN histórico no se fuerza: incluye `SourceRegressorIndex=673, IsQ=true`, que representa la Q de `abs(x)` y se anula sin PN. Mantener su coeficiente y declarar 310 grados ajustables sería incorrecto. No se sustituye esa columna por otra ni se rebautiza un modelo de menor dimensión como P310.

### FLOPs

Se reutilizan `countModelOperations`, `countModelFLOPs` y `getFLOPConvention`: multiplicación/suma real 1/1, compleja 6/2. Los históricos no cambian.

Para IQ se adapta únicamente el calendario analítico: se conserva el coste compartido de magnitudes y potencias, se forma solo cada componente real seleccionada de los términos no lineales, y se contabilizan `2F` multiplicaciones de coeficientes y `2(F-1)` sumas de acumulación. No hay rotación ni restauración de fase. No se asignan los FLOPs de PN-IQ por coincidencia de P.

GMP emplea 155 coeficientes complejos: solo su ponderación/acumulación cuesta `4*155+2*155+2*154=1238` FLOPs. IQ y PN-IQ con F=155 emplean 618 FLOPs para las dos salidas reales. La generación de características y, en PN, la rotación/restauración se añaden aparte.

El contador histórico PN incluye generación/reutilización de magnitudes, dos divisiones para obtener r, rotaciones de portadoras no actuales y auxiliares conjugados, y restauración final de fase. Las sumas/multiplicaciones sí entran en el núcleo; raíces, divisiones y comparaciones quedan aparte, sin peso equivalente implícito. Extracciones real/imaginaria, conjugación, indexación, memoria y centrado por registro no son FLOPs del núcleo. Es una estimación analítica con reutilización, **no** un perfil del código MATLAB ni latencia de hardware.

## D. Resultados reales P310

Ejecución MATLAB R2026a, 12/10/2026. La primera ejecución completa tardó 292 s; las dos selecciones IQ y sus ajustes habían terminado a los 208 s. No se repitió el barrido ni se entrenó una red neuronal.

| Modelo | Parámetros reales | NMSE identificación (dB) | NMSE complemento (dB) | NMSE completo (dB) | FLOPs/muestra |
|---|---:|---:|---:|---:|---:|
| GMP-DOMP histórico | 310 | -38.610857 | -38.403817 | -38.424040 | 1667 |
| PN-IQ-GMP histórico | 310 | -39.054649 | -38.864418 | -38.883032 | 1023 |
| IQ nativo, DOMP propio | 310 | -38.356661 | -38.042672 | -38.073006 | 863 |
| IQ común, DOMP propio | 310 | -38.211173 | -37.935406 | -37.962152 | 870 |
| PN-IQ, soporte del IQ común | 310 | -38.624760 | -38.383511 | -38.406993 | 1017 |

Mejora en el complemento, definida como `NMSE_referencia - NMSE_modelo`; positivo significa mejor:

| Modelo | Frente a GMP (dB) | Frente a PN-IQ histórico (dB) |
|---|---:|---:|
| GMP histórico | 0 | -0.460601 |
| PN-IQ histórico | +0.460601 | 0 |
| IQ nativo | -0.361145 | -0.821746 |
| IQ común | -0.468411 | -0.929012 |
| PN-IQ con soporte IQ común | -0.020306 | -0.480907 |

![NMSE del complemento para los cinco modelos P310](../results/pniq_ablation/P310_20261012_223405_009/nmse_comparison.png)

Desglose del coste estimado, manteniendo las convenciones históricas:

| Modelo | Regresores fuente únicos | Generación: mult. + sumas | Ponderación/acumulación | Rotación/restauración: mult. + sumas | Raíces aparte | Divisiones aparte |
|---|---:|---:|---:|---:|---:|---:|
| GMP histórico | 155 | 415 + 14 | 1238 | 0 | 14 | 0 |
| PN-IQ histórico | 129 | 319 + 14 | 618 | 48 + 24 | 14 | 2 |
| IQ nativo | 115 | 233 + 12 | 618 | 0 | 12 | 0 |
| IQ común | 124 | 240 + 12 | 618 | 0 | 12 | 0 |
| PN-IQ, soporte IQ común | 124 | 315 + 12 | 618 | 48 + 24 | 12 | 2 |

PN-IQ ahorra 644 FLOPs frente al GMP histórico: 620 corresponden al cambio de ponderación/acumulación y los 24 restantes al balance entre generación y fase. **No significa que evalúe la misma función con 644 operaciones menos**: se comparan bases y soportes diferentes. Además, el GMP histórico usa productos complejos literales mientras PN-IQ emplea el calendario analítico reutilizable del proyecto. IQ nativo es aún más barato, pero su NMSE es peor. Menos FLOPs no implica automáticamente mejor modelo ni menor latencia medida.

## E. Controles de soporte e interpretación de diferencias

- **Mismo soporte, distinta representación:** pasar de IQ común a PN-IQ sobre sus mismos 155 identificadores mejora **0.448105 dB** en el complemento, **0.413587 dB** en identificación y **0.444841 dB** en señal completa. El soporte favorecía a IQ, no a PN. Es evidencia concreta de una ventaja de representación PN-IQ en este control, no de una ventaja de rotar un LS complejo.
- **Misma representación IQ, cambio de candidatos y reselección:** el diccionario nativo mejora **0.107266 dB** respecto al común. Sus soportes solo comparten **75/155** identificadores; la eliminación de candidatos afecta también a la trayectoria de selección. No se interpreta como el efecto aislado de una única columna.
- **PN-IQ con soporte propio frente al soporte IQ:** el histórico mejora **0.480907 dB**. Solo comparte **73/155** identificadores con el IQ común. Aquí cambian el soporte y el candidato PN exclusivo; no se puede llamar a toda esa diferencia «ganancia DOMP» independiente del diccionario.
- **I/Q libre sin fase frente a GMP:** IQ nativo empeora **0.361145 dB** con el mismo P. No hay evidencia de que liberar I/Q por sí solo explique la superioridad de PN-IQ; el presupuesto igual no hace anidadas las clases dispersas.

Las diferencias forman una descomposición aritmética dependiente de la ruta elegida, no porcentajes universales de contribución de «fase», «I/Q» y «DOMP». No se ha ajustado ninguna variante a partir del NMSE del complemento.

### Rango y condicionamiento

Los tres ajustes nuevos tienen rango normalizado **155/155**, de modo que ajustan 310 coeficientes reales sin columnas nulas. Los números de condición 2, tras normalizar por pico, son:

| Ajuste | Condición 2 | Repetición LS, error relativo de predicción de identificación |
|---|---:|---:|
| IQ nativo | `1.054553e8` | 0, igualdad exacta |
| IQ común | `1.128645e8` | 0, igualdad exacta |
| PN-IQ en soporte IQ común | `1.176901e8` | No se repitió el ajuste; sí su inferencia guardada |

Son matrices de rango completo pero mal condicionadas. El ajuste principal sigue siendo LS sin Ridge por diseño, no porque se haya demostrado óptimo para robustez o cuantización. Los resultados son reproducibles en doble precisión; no garantizan robustez en precisión simple/fija. Las cancelaciones entre coeficientes impiden interpretar cada peso como una contribución física independiente.

## F. Equivariancia y términos auxiliares

Se rotó la **captura original** a 0°, 45°, 90°, 180° y 270°, se eliminó de nuevo su DC siguiendo ONLINE y se conservaron coeficientes, normalizaciones físicas y taps periódicos. En esta captura no hay muestras de entrada centrada exactamente cero. Se calculó `E(theta)=norm(f(phase*x)-phase*f(x))/norm(f(x))` en doble precisión, sin reajustar nada.

### Predictores completos

| Rotación | GMP histórico | PN-IQ histórico | IQ nativo | IQ común | PN-IQ, soporte IQ común |
|---:|---:|---:|---:|---:|---:|
| 0° | 0 | 0 | 0 | 0 | 0 |
| 45° | `2.17736e-3` | `2.56166e-3` | `4.28693e-3` | `4.68995e-3` | `1.79194e-3` |
| 90° | `3.09042e-3` | `3.65358e-3` | `6.06263e-3` | `6.63260e-3` | `2.53699e-3` |
| 180° | `5.75538e-4` | `9.66362e-4` | `5.88943e-12` | `6.59783e-12` | `9.42625e-12` |
| 270° | `3.09025e-3` | `3.65082e-3` | `6.06263e-3` | `6.63260e-3` | `2.53699e-3` |

### Qué auxiliares se seleccionaron

En la población global, 672 es `conj(x)` y 673 es `abs(x)`:

| Modelo | Fuente / componente | Coeficiente físico complejo o par I/Q |
|---|---|---|
| GMP | 672 | `-0.0017088 + j*0.00016104` |
| GMP | 673 | `-0.00032248 - j*0.0000036279` |
| PN-IQ histórico | 672, Q | `(-0.0011050, -0.0016307)` |
| PN-IQ histórico | 672, I | `(-0.0018901, -0.00079351)` |
| PN-IQ histórico | 673, Q | `(0.00012405, -0.00075477)` |
| IQ nativo | Ninguno explícito | — |
| IQ común | 672, Q | `(-122.99, 242.89)` |
| PN-IQ, soporte IQ común | 672, Q | `(-0.0011775, -0.0016259)` |

La tabla muestra coeficientes en su escala física, no coeficientes de columnas unitarias. Para valorar su influencia se calcula la norma de su contribución a la salida, dividida por la norma de la predicción completa:

| Modelo | Contribución conjugada relativa | Contribución de envolvente relativa |
|---|---:|---:|
| GMP histórico | `0.001531515` | `0.0002877691` |
| PN-IQ histórico | `0.001794153` | `0.0004831809` |
| IQ nativo | 0 | 0 |
| IQ común | **171.9850** | 0 |
| PN-IQ, soporte IQ común | `0.001268493` | 0 |

El valor grande del IQ común es una **cancelación entre términos**, no una salida final 172 veces mayor: `Im(conj(x))=-Im(x)` sustituye a una componente canónica excluida del conjunto común, y su aportación se compensa con las otras columnas. La etiqueta «auxiliar» depende de la base. Eliminar ese término del IQ común destruye esa compensación; no diagnostica una pequeña perturbación aislada ni justifica podarlo.

### Intervención sin reentrenamiento

Se restaron de las predicciones las contribuciones de los dos auxiliares con los coeficientes congelados. El error del predictor intervenido usa su propia predicción a 0° como referencia y denominador. Además, `UnexplainedErrorRelativeToBaseline` guarda la diferencia vectorial que no explican los auxiliares usando el denominador del baseline, evitando confundir un cambio de denominador con una explicación causal.

| Modelo | E(90°) completo | E(90°) sin auxiliares, pesos congelados |
|---|---:|---:|
| GMP histórico | `3.09042e-3` | `8.3896e-12` |
| PN-IQ histórico | `3.65358e-3` | `1.0473e-11` |
| IQ nativo | `6.06263e-3` | `6.06263e-3` |
| IQ común | `6.63260e-3` | **1.407499** |
| PN-IQ, soporte IQ común | `2.53699e-3` | `8.3443e-12` |

![Equivariancia y contribución de auxiliares en los tres modelos principales](../results/pniq_ablation/P310_20261012_223405_009/equivariance_auxiliaries.png)

En GMP y PN-IQ históricos, los auxiliares explican el defecto de equivariancia hasta un residuo de aproximadamente `1e-11`, también con el denominador original. Esta escala es compatible con redondeo amplificado por los coeficientes/condicionamiento, y está muy por debajo del defecto de `1e-3`. A 180° la contribución conjugada cambia con el signo correcto y la envolvente no: se observa exactamente la relación prevista `E(180°)=2*RMS_relativo(envolvente)` dentro del error numérico.

En IQ nativo no había auxiliares explícitos que retirar: la falta de equivariancia procede de la libertad I/Q sobre componentes de términos canónicos, es decir, de su rama conjugada implícita. En PN-IQ esa libertad no rompe por sí misma la equivariancia de los términos canónicos; el defecto medido procede de los auxiliares no canónicos.

No se ha reentrenado un modelo sin auxiliares, ni medido su eventual NMSE óptimo. La intervención es solo causal/algebraica sobre pesos congelados. No se eliminan los términos del baseline: la captura podría contener asimetrías que un modelo estrictamente equivariante no represente, cuestión que esta prueba no resuelve.

## G. Qué se puede atribuir y qué no

**A — ¿PN sola mejora un GMP complejo?** No crea capacidad nueva bajo la transformación común de filas. Se ha demostrado algebraicamente y verificado con LS/DOMP pequeños. No se presenta como arquitectura distinta.

**B — ¿Qué aporta I/Q libre?** Cambia los acoplamientos de coeficientes y el espacio de funciones; no garantiza mejora a P fijo. En este P310, IQ sin PN es peor que GMP, aunque más barato.

**C — ¿Qué aporta combinar PN e I/Q?** Se observa una mejora con selección propia frente a IQ nativo de 0.821746 dB, y una mejora de 0.448105 dB sobre los mismos identificadores del soporte IQ común. Existe evidencia de ventaja de representación en ese control concreto.

**D — ¿Representación, selección o diccionario?** Interactúan. El soporte fijo demuestra que no todo proviene de reseleccionar columnas; la diferencia de 0.480907 dB hasta el PN-IQ histórico demuestra que ese soporte fijo tampoco explica toda su ventaja. No se ha aislado un porcentaje causal para DOMP frente a disponibilidad de candidatos.

**E — ¿Por qué menos FLOPs?** Se seleccionan características reales con dos ponderaciones reales, en lugar de 155 productos coeficiente-regresor complejos. Cambian también los regresores fuente únicos, su generación y el aprovechamiento de identidades PN. Las cifras son calendarios analíticos, no benchmarks de MATLAB o FPGA.

**F — ¿Auxiliares y equivariancia?** La explicación algebraica se confirma cuantitativamente en los baselines. Pero en IQ sin PN también existe una rama conjugada implícita y las etiquetas de los átomos no constituyen una descomposición física única.

## H. Validación, reproducibilidad y límites

### Comprobaciones completadas

- `run_pniq_ablation_test`: PASS. Fixture sintético de 256 muestras y P12, solo modelos lineales. Comprueba equivalencia LS compleja, DOMP bajo rotación común, objetivo complejo frente a dos salidas reales, identidad widely-linear, conteos de candidatos, rango, presupuesto y soporte común, reconstrucción/restauración y recarga.
- Ambos predictores históricos P310: **igualdad exacta** frente a las señales históricas guardadas en el diagnóstico anterior. NMSE recalculadas compatibles con las métricas del bundle con tolerancia `1e-8 dB`; pequeñas diferencias de redondeo entre evaluación por bloques y ajuste sobre identificación no modifican los archivos originales.
- Los cinco modelos guardados: **igualdad exacta de todas sus predicciones al recargarlos en otra sesión MATLAB**. No se repitió DOMP para esta comprobación.
- Repetición de LS sobre los dos soportes IQ ya elegidos: predicciones de identificación idénticas, error relativo 0. No se repitieron las selecciones.
- Los tres soportes nuevos: rango 155, presupuesto 310; en el control recíproco coinciden exactamente los 155 pares `(SourceRegressorIndex, IsQ)`.
- Rotaciones comprobadas con la captura original y DC ONLINE; no se recalibran coeficientes ni escala.
- `checkcode`: sin mensajes en los tres archivos MATLAB nuevos. `git diff --check`: sin errores; comprobados también los archivos nuevos sin seguimiento.
- Dos PNG revisados visualmente. No se generan diez figuras ni se regenera ninguna histórica.

No se ha ejecutado una batería neuronal ni el probe PA. No se presentan como superados tests que no se hayan ejecutado.

### Archivos y comandos

Nuevos archivos, sin cambios en funciones científicas históricas:

- `runPNIQAblation.m`: coordinación, dos selecciones IQ, LS y controles, costes y figuras.
- `toolbox/sweep/predictPNIQAblation.m`: único helper de inferencia experimental por bloques, reutilizado para recarga y diagnóstico; no modifica `applyDPDModels`.
- `tests/run_pniq_ablation_test.m`: prueba pequeña.
- Este informe.

Resultados reales: [carpeta P310](../results/pniq_ablation/P310_20261012_223405_009/). Contiene `ablation.mat` (configuración, split, modelos, coeficientes, soportes, predicciones y tablas), `metrics.csv`, `rotations.csv`, `auxiliary_contributions.csv`, `numerical_checks.csv` y los dos PNG. El MAT ocupa aproximadamente 43 MB; no guarda matrices de regresión. Los resultados sintéticos están separados en carpetas P12.

Para consultar, sin volver a ajustar:

```matlab
S = load('results/pniq_ablation/P310_20261012_223405_009/ablation.mat');
S.result.metrics
S.result.rotations
```

Para repetir voluntariamente el experimento en otra carpeta nueva:

```matlab
result = runPNIQAblation( ...
    'results/dpd_sweep/tp98e6f7ac_f3cd_44c7_9ce3_fd8619319748', 310);
```

Para inferir nuevamente un modelo experimental ya guardado, sin ajuste:

```matlab
addpath(genpath('toolbox'));
D = load(fullfile(S.result.sourceDirectory,'experiment.mat'));
yhat = predictPNIQAblation(S.result.models{3}, D.x, D.cfg.gmp.blockSize);
```

El helper recibe la señal **ya preprocesada**; no es una nueva interfaz RF. Los bundles ONLINE P310/P340 y todo el WIP previo permanecen intactos. No se hizo commit ni push.

### Limitaciones

Una única captura ILC. El estudio inicial de P310 se amplía a P100 y P500 en la sección J. El complemento procede de la misma captura y P310 fue seleccionado usándolo; no es un test independiente. No hay nuevas medidas RF, ACLR, ACPR o EVM. Los candidatos nativos difieren y el soporte fijado favorece a IQ; falta una comparación neutral sobre un soporte predefinido o selección en otra captura. Las matrices están mal condicionadas y no se ha estudiado precisión reducida. Los FLOPs son estimaciones, no latencia ni consumo.

La ampliación P100/P500 estudia dependencia de P, no constituye una confirmación independiente entre capturas. No se repite P310 ni se lanzan los 49 presupuestos.

## I. Conclusión y siguiente experimento

La superioridad observada de PN-IQ no puede atribuirse a «rotar GMP»: ese control complejo es equivalente. Tampoco la explica por sí sola liberar I/Q, que aquí empeora el resultado con P fijo. Hay una mejora legítima y medible de **0.448105 dB con soporte idéntico** al combinar fase e I/Q, además de efectos relevantes de selección y candidatos que no se pueden separar completamente con estos controles.

Los auxiliares explican la pérdida de equivariancia de los baselines; no se ha demostrado que eliminarlos mejore NMSE ni linealización RF. IQ sin PN presenta además una pérdida de equivariancia intrínseca a su rama conjugada libre.

El siguiente experimento de mayor valor es repetir un protocolo P310 predefinido en **otra captura ILC independiente**, manteniendo los controles de soporte y separando selección/ajuste de evaluación externa. Puede incluirse un soporte predefinido en la primera captura para evaluar cuánto depende la ventaja de DOMP y de la captura. La dependencia del presupuesto dentro de esta captura se estudia a continuación; el rendimiento RF requiere autorización y medidas del laboratorio. No se justifica aún una afirmación universal ni repetir automáticamente todo el barrido.

## J. Cierre experimental P100 / P310 / P500 (2026-10-13)

### J1. Ejecución y resultados conjuntos

P100 terminó en 37.0 s y P500 en 485.3 s. Se reutilizaron los bundles históricos, la captura DC-removed, las 49151 filas de identificación, las 442369 de evaluación y los taps periódicos. P100 usa los primeros 50 identificadores de las dos rutas IQ guardadas en P310. P500 selecciona 250; sus primeros 155 identificadores coinciden **exactamente** con P310 en ambas rutas. No se repitió DOMP histórico ni el experimento P310, no se entrenó PNNN y no se utilizó Ridge.

Resultados nuevos: [P100](../results/pniq_ablation/P100_20261012_232926_441/) y [P500](../results/pniq_ablation/P500_20261012_233006_611/). Cada carpeta conserva modelos, predicciones, métricas y los diagnósticos de fase que ya produce la función. No se añadieron figuras: las tablas permiten comparar los tres presupuestos sin rehacer las figuras P310.

| P | Modelo | NMSE identificación (dB) | NMSE evaluación (dB) | NMSE completa (dB) | FLOPs/muestra |
|---:|---|---:|---:|---:|---:|
| 100 | GMP histórico | -37.714283 | -37.600338 | -37.611576 | 549 |
| 100 | PN-IQ histórico | -37.343365 | -37.262215 | -37.270246 | 364 |
| 100 | IQ nativo | -36.000305 | -35.846605 | -35.861701 | 287 |
| 100 | IQ común | -35.636785 | -35.505158 | -35.518116 | 294 |
| 100 | PN-IQ, soporte IQ común | -36.126124 | -36.000361 | -36.012749 | 378 |
| 310 | GMP histórico | -38.610857 | -38.403817 | -38.424040 | 1667 |
| 310 | PN-IQ histórico | -39.054649 | -38.864418 | -38.883032 | 1023 |
| 310 | IQ nativo | -38.356661 | -38.042672 | -38.073006 | 863 |
| 310 | IQ común | -38.211173 | -37.935406 | -37.962152 | 870 |
| 310 | PN-IQ, soporte IQ común | -38.624760 | -38.383511 | -38.406993 | 1017 |
| 500 | GMP histórico | -38.825805 | -38.539738 | -38.567454 | 2623 |
| 500 | PN-IQ histórico | -39.354424 | -39.075585 | -39.102620 | 1549 |
| 500 | IQ nativo | -38.678558 | -38.277772 | -38.316147 | 1359 |
| 500 | IQ común | -38.614896 | -38.241197 | -38.277077 | 1363 |
| 500 | PN-IQ, soporte IQ común | -39.133556 | -38.838161 | -38.866753 | 1543 |

La tabla conjunta copia las métricas históricas originales después de verificar las recalculadas con tolerancia de `1e-8 dB`. Así conserva exactamente los números del bundle, sin confundir redondeos de evaluación por bloques con cambios de modelo. El MAT P310 original no se reescribe.

Mejora positiva = menor NMSE de evaluación con PN; diferencia FLOPs = coste PN menos coste del comparador:

| P | PN histórico vs GMP (dB) | PN histórico vs IQ nativo (dB) | PN vs IQ, soporte común (dB) | ΔFLOPs PN−GMP | ΔFLOPs PN−IQ nativo | ΔFLOPs común PN−IQ |
|---:|---:|---:|---:|---:|---:|---:|
| 100 | -0.338123 | 1.415610 | 0.495203 | -185 | +77 | +84 |
| 310 | 0.460601 | 0.821746 | 0.448105 | -644 | +160 | +147 |
| 500 | 0.535847 | 0.797812 | 0.596964 | -1074 | +190 | +180 |

La ventaja de representación del control con los mismos identificadores **se mantiene en los tres presupuestos**, pero no crece monótonamente. En cambio, la superioridad del PN-IQ histórico sobre GMP **depende del presupuesto**: GMP gana en P100 y PN-IQ en P310/P500. PN-IQ histórico supera a IQ nativo en los tres casos, con coste adicional; es siempre más barato que GMP. Esto no atribuye toda la ganancia a una rotación: intervienen la libertad de coeficientes I/Q, el espacio de funciones, los candidatos y la selección, con interacciones entre ellos. El soporte común se seleccionó para IQ, no es un soporte neutral.

### J2. Diagnóstico numérico: método

`diagnosePNIQAblation` lee los modelos guardados sin modificarlos. Reconstruye exclusivamente las columnas seleccionadas sobre identificación: complejas para GMP y reales para IQ/PN-IQ, estas últimas con sus mismos pares `(SourceRegressorIndex, IsQ)`. Emplea exactamente el escalado del LS: `B=A./max(abs(A),[],1)` y objetivo de pico unitario, rotado para PN-IQ. No utiliza una matriz de Gram ni `cond` sobre el diccionario completo.

La SVD económica de B da los extremos singulares, `kappa2=sigmaMax/sigmaMin` y el rango con `max(size(B))*eps(sigmaMax)`, el criterio del ajuste existente. El rango complejo de GMP y el rango real de IQ/PN cuentan columnas en sus respectivos cuerpos; cada columna seleccionada lleva dos parámetros reales. Se marca casi dependencia cuando `sigmaMin/sigmaMax<1e-8`; es un indicador explícito, no una prueba automática de predicciones inestables. También se registra la razón entre las normas 2 máxima y mínima de las columnas ya escaladas por pico, para distinguir disparidad residual de escalas de dependencia entre características.

La comprobación dirigida repite LS con soporte congelado y verifica repetibilidad exacta. Después perturba únicamente el objetivo normalizado con un vector de norma relativa `1e-12`, orientado según la dirección singular izquierda más débil. Mide cambio de coeficientes, predicción completa y NMSE de evaluación respecto al refit sin perturbación. Se verifican asimismo las predicciones de los modelos recargados, exactamente iguales a las guardadas. Ninguno de estos ajustes de diagnóstico reemplaza coeficientes ni resultados históricos; no hay regularización nueva ni nuevas selecciones.

Los 15 modelos tienen exactamente P parámetros reales y los controles comunes conservan exactamente la correspondencia de características. No aparece ninguna incompatibilidad estructural en esos controles. Transferir directamente el soporte PN histórico a IQ sería otra cuestión: P100 no contiene características inadmisibles en IQ, pero P310 y P500 contienen una, la componente Q del auxiliar de envolvente, nula sin PN. No se fuerza esa transferencia ni se modifica el soporte para favorecer una comparación.

### J3. Condicionamiento y sensibilidad observados

| P | Modelo | Rango / columnas | Sigma máxima | Sigma mínima | Condición 2 |
|---:|---|---:|---:|---:|---:|
| 100 | GMP histórico | 50 / 50 | 200.120 | 2.16953e-6 | 9.22413e7 |
| 100 | PN-IQ histórico | 50 / 50 | 160.202 | 2.71584e-5 | 5.89879e6 |
| 100 | IQ nativo | 50 / 50 | 114.041 | 6.22674e-3 | 1.83147e4 |
| 100 | IQ común | 50 / 50 | 117.258 | 4.47133e-3 | 2.62244e4 |
| 100 | PN-IQ, soporte IQ común | 50 / 50 | 148.948 | 5.14970e-3 | 2.89236e4 |
| 310 | GMP histórico | 155 / 155 | 257.924 | 7.81626e-7 | 3.29984e8 |
| 310 | PN-IQ histórico | 155 / 155 | 209.070 | 1.57178e-6 | 1.33015e8 |
| 310 | IQ nativo | 155 / 155 | 159.986 | 1.51710e-6 | 1.05455e8 |
| 310 | IQ común | 155 / 155 | 171.274 | 1.51752e-6 | 1.12864e8 |
| 310 | PN-IQ, soporte IQ común | 155 / 155 | 185.014 | 1.57204e-6 | 1.17690e8 |
| 500 | GMP histórico | 250 / 250 | 279.358 | 7.78927e-7 | 3.58645e8 |
| 500 | PN-IQ histórico | 250 / 250 | 228.216 | 1.56972e-6 | 1.45386e8 |
| 500 | IQ nativo | 250 / 250 | 175.825 | 1.51572e-6 | 1.16001e8 |
| 500 | IQ común | 250 / 250 | 180.814 | 1.51555e-6 | 1.19306e8 |
| 500 | PN-IQ, soporte IQ común | 250 / 250 | 213.054 | 1.57013e-6 | 1.35692e8 |

Todos tienen rango numérico completo. El indicador de casi dependencia se activa para los diez modelos P310/P500; ninguno de P100 cruza ese umbral, aunque GMP P100 ya está cerca y no debe calificarse de bien condicionado. La razón de normas 2 de columnas tras el escalado de pico está entre 9.83 y 15.58: la disparidad residual de escalas no explica por sí sola condiciones de hasta `3.59e8`.

**PN-IQ no está mejor condicionado que IQ en estos controles.** Con el mismo soporte, su condición es aproximadamente 10.3%, 4.3% y 13.7% mayor en P100/P310/P500. El PN-IQ histórico también tiene peor condición que IQ nativo, aunque mejor que GMP histórico. El deterioro con P es especialmente acusado en IQ entre P100 y P310; de P310 a P500 es más moderado en todas las variantes. Por tanto, la ventaja de NMSE de PN-IQ no puede explicarse aquí mediante un mejor número de condición.

Máximos absolutos entre los cinco modelos de cada presupuesto en la perturbación dirigida del objetivo:

| P | Cambio relativo de coeficientes | Cambio relativo de predicción completa | Cambio absoluto NMSE evaluación (dB) |
|---:|---:|---:|---:|
| 100 | 4.314e-10 | 1.330e-11 | 9.372e-12 |
| 310 | 1.067e-9 | 3.019e-11 | 1.599e-10 |
| 500 | 1.059e-9 | 3.998e-11 | 1.835e-10 |

Los refits sin perturbación reproducen exactamente las predicciones guardadas en los 15 casos; también lo hace la recarga de cada modelo. Repetir LS con las mismas entradas produce exactamente los mismos coeficientes. No hay evidencia en estas pruebas de que las diferencias de NMSE, incluidas las pequeñas de centésimas de dB, se deban a inestabilidad numérica observable. Las mejoras de 0.448–0.597 dB del control común están muchos órdenes por encima de los efectos medidos.

Esto **no** demuestra robustez frente a cualquier perturbación: solo se perturbó el objetivo, no la matriz, y se utilizó doble precisión en esta captura. No se estudió precisión reducida, ruido de medida adicional, cambios de captura ni confianza estadística. Las cancelaciones entre coeficientes descritas en P310 siguen desaconsejando interpretar sus magnitudes como contribuciones físicas independientes. El condicionamiento merece prudencia para implementación y cuantización, pero no invalida por sí solo las comparaciones observadas.

### J4. Artefactos, validación y cierre

La [carpeta conjunta](../results/pniq_ablation/summary_20261013_164459/) contiene `metrics.csv` (15 filas con los tres dominios NMSE y FLOPs), `advantages.csv`, `conditioning.csv` (incluye escalas y sensibilidad completas) y `comparison.mat`. Los tres experimentos originales permanecen en sus carpetas, sin sobrescritura. Para repetir únicamente una auditoría, no el entrenamiento:

```matlab
addpath(genpath('toolbox'));
d = diagnosePNIQAblation( ...
    'results/pniq_ablation/P500_20261012_233006_611/ablation.mat');
```

Cambios de esta ampliación: `runPNIQAblation.m` acepta prefijos IQ guardados y permite omitir figuras individuales; conserva exactamente las métricas históricas verificadas. Se amplió `tests/run_pniq_ablation_test.m` para comprobar la reutilización y esas métricas. Se añadió únicamente `toolbox/sweep/diagnosePNIQAblation.m` para la auditoría reutilizable descrita arriba y se amplió este informe. `predictPNIQAblation.m` y las funciones históricas de entrenamiento/inferencia no se modificaron.

Validación completada: test sintético existente PASS (incluida reutilización de rutas), prefijos reales P310/P500 idénticos, presupuestos y mapas comunes correctos, baselines conservados, recarga/refit/repetición LS exactos en los 15 modelos, `checkcode` sin mensajes y `git diff --check` sin errores. Se reutilizó el test ya ejecutado antes de la interrupción; no se repitieron cálculos terminados. No se entrenó PNNN, no se ejecutó el barrido completo, no hubo commit ni push.

**Conclusión:** se cierra la ablación en tres presupuestos con evidencia consistente a favor de la representación PN-IQ sobre el soporte común de esta captura. La superioridad sobre GMP histórico depende de P y aparece en P310/P500, no en P100. Las pruebas dirigidas no atribuyen esas ventajas a artefactos numéricos, pero tampoco las convierten en una afirmación universal o en una ganancia exclusivamente causal de la fase. El siguiente paso recomendado sigue siendo validar un protocolo P310 predefinido en otra captura ILC independiente, separando selección y evaluación externa; no otro barrido completo de esta captura.
