# Equivariancia de fase y diccionario GMP restringido

Fecha: 20/10/2026. Ejecutor `DPD/experiments/runGMPPhaseStudy.m`.
Directorio separado: `results/gmp_phase_study/20261020_205246/`.
Captura ILC y split: los mismos de
`results/dpd_sweep/tp98e6f7ac_f3cd_44c7_9ce3_fd8619319748/`.
No cambia el baseline, P310 ni ninguna PNNN.

## Propiedad matemática verificada en la implementación

`buildGMPRegressorRows` evalúa productos de taps de x, taps conjugados y
envolventes. Para cada monomio, incluidos retardos periódicos,
q = numel(X) - numel(Xconj). Entonces
g(exp(j phi)x) = exp(j q phi)g(x). Las envolventes no aportan fase.
La clasificación se obtiene de esos factores, no del nombre del regressor.

| Familia real | Candidatos | q | Equivariancia de grado 1 |
|---|---:|---:|---|
| x actual | 1 | 1 | Sí |
| x retardada | 10 | 1 | Sí |
| GMP portadora × productos de envolventes | 660 | 1 | Sí |
| conj(x) actual | 1 | -1 | No |
| abs(x) actual | 1 | 0 | No |

`regressorManager.finalizePopulation` añade x, conj(x), abs(x) y deduplica:
x ya existía; los dos últimos quedan como auxiliares 672 y 673.
El diccionario restringido tiene 671 columnas complejas y 1281 candidatos
PN-IQ reales, frente a 673/1285. Se conservan los 61 ceros Q estructurales
excluidos, el orden relativo y todas las fórmulas de los candidatos restantes.

Para PN-IQ, escribiendo r=conj(x)/abs(x) y combinando los coeficientes I/Q,
la salida tiene la forma Uc + conj(r)^2 conj(U)d. Sus grados de fase son
q y 2-q: los canónicos q=1 son equivariantes aun con coeficientes I/Q libres;
conj(x) introduce grados -1 y 3, y abs(x) grados 0 y 2.
Rotar conjuntamente matriz y objetivo de LS preserva su norma residual,
pero NO prueba esta propiedad al reconstruir características desde otra entrada.

La afirmación PN-IQ presupone x actual no nula: en x(n)=0 el código fija r=1.
Con taps retardados no nulos y coeficientes I/Q libres, esa convención puede
romper la equivariancia incluso sin auxiliares. Por ejemplo, para una sola
característica Re(r*x(n-1)), x(n)=0 y x(n-1)=1, rotar 90° da 0 en vez de j.
No se modifica la convención histórica. En muestras pequeñas pero no nulas
la identidad exacta sigue siendo válida, con sensibilidad de coma flotante.

La simetría es una hipótesis del modelo de banda base, no una ley garantizada
para cualquier cadena de adquisición/ILC. Desbalances I/Q u otros mecanismos
pueden motivar una respuesta no equivariante, pero estos experimentos no los
identifican. «Auxiliar no equivariante» no significa «término físicamente
erróneo»; por eso se mide también la pérdida o ganancia de NMSE al excluirlo.

## Diseño controlado y recuperable

Primero se diagnostican los modelos históricos congelados P100/P310/P500.
Se reutilizan sus predicciones y las rotaciones 0/90/180° de la ablación
anterior; únicamente se calculan 30/60° que faltaban. El campo `Reused`
identifica la procedencia de cada fila. Se verifica igualdad de modelos
con los bundles del sweep, sin volver a ajustarlos.

Después se calculan dos trayectorias DOMP NUEVAS sobre los candidatos
restringidos, de 250 átomos cada una. La selección es secuencial: el prefijo
de 155 no depende del número de iteraciones posteriores. Se ajusta y valida
P310 primero; después P100 y P500 reutilizan esos prefijos. No se usa ninguna
trayectoria del diccionario completo. Se conservan pico por columna, objetivo
unit-peak, LS `lsqminnorm`, tolerancia de rango, muestras, bloques y conteo.

`restricted_path_1/2.mat` y `restricted_P0100/0310/0500.mat` son checkpoints
solo de este experimento. Cada punto guarda modelos, predicciones, soportes
con identificadores del diccionario original, métricas y rotaciones.
La reanudación explícita es `runGMPPhaseStudy('ruta_de_esta_ejecucion')`.
No se escriben archivos en el sweep ni en la ablación histórica.

El diagnóstico rota la captura raw, vuelve a quitar DC como ONLINE y usa
norm(f(rotated)-phase*f(original))/norm(f(original)) en el registro completo.
La evaluación NMSE usa únicamente el complemento histórico. Se distinguen
la supresión congelada de auxiliares y el reentrenamiento lineal restringido:
son intervenciones diferentes. Las contribuciones relativas no son porcentajes
causales ni energías aditivas, porque existen términos cruzados.

## Resultados congelados verificados

| P | GMP: fuentes conj / env | PN-IQ: características conj / env | E(90°) GMP | E(90°) PN-IQ |
|---:|---:|---:|---:|---:|
| 100 | 1 / 0 | 2 / 0 | 0.00307587 | 0.00357259 |
| 310 | 1 / 1 | 2 / 1 | 0.00309042 | 0.00365358 |
| 500 | 1 / 1 | 2 / 1 | 0.00310159 | 0.00367620 |

El CSV de composición contiene además x actual, retardadas y no lineales,
contando fuentes únicas y características por separado. PN-IQ puede usar I y Q
de una misma fuente; no se debe sumar esas dos columnas como dos monomios
complejos independientes.

En P310 las normas relativas de contribución conjugada/envolvente son
0.00153151/0.000287769 para GMP y 0.00179415/0.000483181 para PN-IQ.
No representan energías sumables ni la caída de NMSE al retirar cada término.
PN-IQ histórico tiene MENOR error de imitación pero MAYOR defecto de
equivariancia que GMP en las rotaciones examinadas: no es más equivariante
por llamarse phase-normalized.

En P310, E(30°)/E(60°)/E(180°) es
0.00153895/0.00266860/0.000575538 para GMP y
0.00180981/0.00314264/0.000966362 para PN-IQ.
En P100 E(180°) cae a ~1e-11/1e-13 porque no se seleccionó la envolvente:
examinar solo 180° habría ocultado el defecto debido al conjugado.

La supresión congelada de ambos auxiliares ya guardada en la ablación
reduce E(90°) a 8.23e-12/6.89e-14 (P100),
8.39e-12/1.05e-11 (P310) y 8.44e-12/1.05e-11 (P500).
Esto localiza el defecto de fase, pero no predice el NMSE de una nueva
selección restringida. Los 12 valores nuevos de 30/60° se calcularon una
sola vez; los 18 restantes se reutilizaron con su procedencia explícita.

## Diccionario restringido: resultados completados

Los seis ajustes terminaron el 20/10/2026 a las 21:26. Al reanudar el
21/10 se recuperaron `study.mat` y los tres checkpoints: **no se repitieron
DOMP, LS, SVD ni inferencias**. El proceso original seguía vivo sin producir
archivos desde la exportación de las 21:27; se cerró exclusivamente ese PID
tras identificar su comando. `complete_exports.m`, conservado junto a los
resultados, completa solo las exportaciones ausentes desde los datos guardados.

NMSE sobre las mismas 442369 muestras del complemento de evaluación;
Δ positivo significa empeoramiento respecto al histórico:

| P | Modelo | Histórico (dB) | Restringido (dB) | Δ (dB) | FLOPs histórico → restringido |
|---:|---|---:|---:|---:|---:|
| 100 | GMP | -37.600338 | -37.560480 | +0.039858 | 549 → 554 |
| 100 | PN-IQ | -37.262215 | -37.295086 | -0.032871 | 364 → 358 |
| 310 | GMP | -38.403817 | -38.329724 | +0.074093 | 1667 → 1663 |
| 310 | PN-IQ | -38.864418 | -38.714206 | +0.150212 | 1023 → 1026 |
| 500 | GMP | -38.539738 | -38.430443 | +0.109295 | 2623 → 2628 |
| 500 | PN-IQ | -39.075585 | -38.970250 | +0.105335 | 1549 → 1591 |

`metrics.csv` conserva además identificación y señal completa, sin sustituir
ninguna métrica del sweep. Los seis modelos tienen exactamente P parámetros
reales activos y P/2 columnas seleccionadas. No se añadió Ridge.

| P | Modelo | Columnas comunes con el histórico / total | Condición 2 histórica → restringida | Rango restringido |
|---:|---|---:|---:|---:|
| 100 | GMP | 49 / 50 | 9.22413e7 → 9.22498e7 | 50 |
| 100 | PN-IQ | 48 / 50 | 5.89879e6 → 1.01236e8 | 50 |
| 310 | GMP | 104 / 155 | 3.29984e8 → 1.18100e8 | 155 |
| 310 | PN-IQ | 94 / 155 | 1.33015e8 → 1.32469e8 | 155 |
| 500 | GMP | 168 / 250 | 3.58645e8 → 3.89285e8 | 250 |
| 500 | PN-IQ | 151 / 250 | 1.45386e8 → 2.58934e8 | 250 |

El soporte se compara por identificador original y componente I/Q, no por
posición local del diccionario. `supports.csv` conserva cada selección y orden;
`classification.csv` permite recuperar su familia. Excluir solo dos fuentes
puede alterar muchas decisiones DOMP posteriores: en P310 cambia 51/155
columnas GMP y 61/155 PN-IQ. No es equivalente a anular dos coeficientes.
La condición se calcula sobre columnas escaladas por pico, como en el ajuste,
con SVD económica y tolerancia `max(size(A))*eps(sigmaMax)`. Las condiciones
históricas se reutilizan del diagnóstico anterior, no se recalculan.

| P | Máximo E restringido GMP | Máximo E restringido PN-IQ |
|---:|---:|---:|
| 100 | 1.48035e-11 | 7.85615e-12 |
| 310 | 1.49148e-11 | 1.34795e-11 |
| 500 | 1.49848e-11 | 1.34184e-11 |

Máximos sobre 0/30/60/90/180 grados. La recarga/reproducción de cada modelo
coincide exactamente con su predicción de ajuste (`ReplayRelativeError=0`).
La entrada centrada no contiene ceros exactos; su módulo mínimo es
5.74259118704624e-5. Por tanto, la excepción matemática PN-IQ en x(n)=0 no se
activa aquí. Las diferencias de fase ~1e-11 son compatibles con redondeo y
cancelación en doble precisión, con matrices de condición ~1e8; no se exige
igualdad bit a bit bajo una transformación aritmética diferente.

## Interpretación y límites

1. **La restricción restaura la equivariancia numérica examinada**, reduciendo
   defectos históricos del orden de 0.003 a ~1e-11. Lo explican los grados de
   fase de los auxiliares, no el nombre de la arquitectura.
2. **No mejora sistemáticamente el NMSE.** Empeoran cinco de seis ajustes
   (0.040–0.150 dB); PN-IQ P100 mejora 0.033 dB. La propiedad algebraica y la
   precisión de imitación son objetivos distintos. Tampoco mejora siempre el
   condicionamiento: PN-IQ P100 pasa de 5.90e6 a 1.01e8.
3. **La ventaja PN-IQ a igual P permanece donde ya existía:** en el diccionario
   restringido supera a GMP en 0.384481 dB (P310) y 0.539807 dB (P500), con
   38.30% y 39.46% menos FLOPs respectivamente. En P100 pierde 0.265394 dB,
   aunque cuesta menos. Los auxiliares no explican por sí solos esa ventaja.
4. El resultado nuevo es una cuantificación controlada del compromiso
   simetría/precisión y del cambio de soportes **en esta captura**, no un
   teorema novedoso ni una innovación DPD demostrada en RF. No justifica
   reemplazar el baseline histórico ni su selección P310.
5. No se ha identificado la causa física de los auxiliares útiles, ni se han
   estimado intervalos de confianza entre capturas. La evaluación pertenece
   a la misma captura y protocolo de desarrollo; faltan otra captura y
   mediciones RF de linealización. Las cinco rotaciones no sustituyen la
   demostración algebraica ni cubren errores de precisión reducida.

Figuras basadas exclusivamente en los resultados guardados:

![NMSE histórico y restringido](../results/gmp_phase_study/20261020_205246/restricted_nmse.png)

![Equivariancia P310](../results/gmp_phase_study/20261020_205246/phase_equivariance_P310.png)

Cada figura dispone de FIG editable, PNG y TikZ. La recuperación no sobrescribe
los FIG/PNG de NMSE ya terminados ni los modelos/checkpoints. El comando de
reanudación del ejecutor devuelve `study.mat` directamente cuando existe;
para este cierre interrumpido se usa exclusivamente `complete_exports.m`.

## Verificaciones del cierre (21/10/2026)

- Tablas de eficiencia idénticas al sweep y splits E/Q idénticos; conteos de
  evaluación, seis presupuestos activos y rangos P/2 correctos; todas las
  predicciones restringidas tienen 491520 muestras.
- Selectores conservados: DPD P310 y PA legacy P340. Los cinco FIG de E/Q
  abren; las dos figuras Q se inspeccionaron visualmente.
- La primera verificación pasó esas aserciones y se detuvo en tres comentarios
  de supresión AGROW obsoletos. Se eliminaron únicamente esos comentarios del
  ejecutor actual; la copia original de ejecución permanece intacta.
- `checkcode` sin mensajes: ambos ejecutores experimentales y los dos scripts
  cortos de recuperación/consulta. Evidencia en `closure_matlab.log`,
  `regularization_matlab.log` y `final_checkcode.log`. El script de la primera
  verificación se conserva como `validate_completed_stages.m`; no se volvió
  a ejecutar lo que ya había pasado.
- 511 archivos históricos con rutas, tamaños y fechas inalterados respecto al
  inventario previo. No se repitieron las pruebas históricas ya superadas.
- La consulta R2026a sin entrenamiento confirmó L2 densa=1e-4, factor pesos=1
  y sesgos=0. El fine-tuning personalizado no añade L2. Se guarda en
  `regularization.mat` y se explica en el plan actualizado.
- No hay asociaciones RF en las 147 filas de `rf_measurements.csv`.
  Las métricas físicas siguen pendientes; no hubo transmisión ni entrenamiento
  neuronal, commit o push.
