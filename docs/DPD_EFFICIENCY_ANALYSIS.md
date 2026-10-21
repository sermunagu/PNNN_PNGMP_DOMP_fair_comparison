# Eficiencia DPD: NMSE, FLOPs y amplitud

Fecha: 20/10/2026. Análisis ejecutado, sin identificación ni entrenamiento nuevo.
Ejecutor: `DPD/experiments/analyzeDPDEfficiency.m`.
Resultados: `results/dpd_efficiency/20261020_204728/`.
Fuente: `results/dpd_sweep/tp98e6f7ac_f3cd_44c7_9ce3_fd8619319748/`.
Captura: `measurements/experiment20260429T134032_xy.mat`, x=u, y=x_ILC.

## Protocolo y trazabilidad

Se reutilizan los 147 puntos (tres familias, P=20:10:500), configuración,
split y predicciones guardadas. No se ejecutó ONLINE. El análisis de amplitud
usa P310, 491520 muestras totales y el mismo complemento de 442369 muestras.
Las predicciones lineales conservan precisión double y PNNN single; las
estadísticas de energía se acumulan en double. La NMSE recalculada con el
helper histórico coincide con la guardada con tolerancia de 1e-8 dB.
Las medias DC son las históricas, no se estiman sobre un subconjunto nuevo.

La evaluación no se usó para ajustar coeficientes finales, pero sí para elegir
P310. Es la misma captura y el split fue informado por el pico del objetivo;
no es un test externo. No se estiman intervalos de confianza tratando las
muestras temporalmente correlacionadas como réplicas independientes.

## Pareto y comparación discreta

Se domina un punto cuando otro tiene FLOPs y NMSE no mayores y alguna
desigualdad estricta. `pareto.csv` contiene indicadores por familia y global.
Las fronteras de familia tienen 49 puntos GMP, 49 PN-IQ y 38 PNNN.
La frontera global tiene **49 PN-IQ, 3 GMP (P20/P30/P40) y ningún PNNN**.
Los otros 95 puntos son dominados globalmente. Esto depende de la captura,
el rango de P y la convención aritmética, no demuestra superioridad universal.

`equal_flops.csv` compara los 145 techos observados, con 435 filas.
Se elige el modelo calculado de menor NMSE que satisface coste <= techo.
No se interpolan modelos. Un tramo horizontal después del último coste de
una familia significa reutilizar su mejor modelo factible, NO extrapolar
un modelo de mayor presupuesto. NaN significa que no hay modelo factible.

| Techo FLOPs | GMP: P / NMSE dB | PN-IQ: P / NMSE dB | PNNN: P / NMSE dB |
|---:|---:|---:|---:|
| 222 | 40 / -34.9197 | 50 / -34.8447 | 30 / -5.0158 |
| 245 | 40 / -34.9197 | 60 / -35.9749 | 40 / -12.2448 |
| 780 | 140 / -37.9939 | 220 / -38.6397 | 310 / -35.4823 |
| 1023 | 180 / -38.1772 | 310 / -38.8644 | 430 / -35.6531 |
| 1667 | 310 / -38.4038 | 500 / -39.0756 | 500 / -35.6946 |

PN-IQ no gana en todos los techos: GMP gana a 222 FLOPs. En la región de
mayor calidad PN-IQ ofrece una ventaja clara: a 780 FLOPs mejora 0.6458 dB
frente a GMP, a 1023 mejora 0.6872 dB y a 1667 mejora 0.6718 dB.
La ventaja a igual coste no debe confundirse con comparar P310 entre familias.
Se comprobó toda la tabla: PN-IQ gana en cada techo observado entre 245 y
2623 FLOPs. Fuera del rango calculado no se hacen afirmaciones. A partir de
1549, su mejor modelo factible es P500; aumentar el techo no crea un modelo nuevo.

En P310: GMP -38.403817 dB/1667 FLOPs; PN-IQ -38.864418 dB/1023;
PNNN -35.482326 dB/780. PN-IQ mejora GMP 0.460601 dB con 38.63% menos
FLOPs. PNNN cuesta menos en ese P pero pierde 3.382092 dB frente a PN-IQ.
En P100, la ablación previa muestra que GMP gana a PN-IQ a igual P:
no se oculta ese resultado.

## Igual calidad objetivo

Modelo realmente calculado de mínimo coste; no son cruces interpolados.

| NMSE objetivo | GMP: P / FLOPs | PN-IQ: P / FLOPs | PNNN: P / FLOPs |
|---:|---:|---:|---:|
| -35 dB | 50 / 275 | 60 / 245 | 240 / 640 |
| -37 dB | 80 / 434 | 100 / 364 | No alcanzado |
| -38 dB | 150 / 812 | 140 / 486 | No alcanzado |
| -38.5 dB | 410 / 2167 | 200 / 686 | No alcanzado |
| -39 dB | No alcanzado | 430 / 1354 | No alcanzado |

Para -38 dB PN-IQ requiere 40.15% menos FLOPs que GMP; para -38.5 dB,
68.34% menos. Son umbrales escogidos para describir los datos, no criterios
de selección nuevos ni promesas de calidad RF.

## Qué cuenta un FLOP

`getFLOPConvention`: suma/producto real=1, producto complejo=6, suma
compleja=2. Los contadores incluyen un calendario de generación de
regresores/características y rotaciones; no son mediciones de instrucciones.

| P310 | Generación | Coeficientes/acumulación | Normalización/restauración de fase | Total |
|---|---:|---:|---:|---:|
| GMP | 429 | 1238 | 0 | 1667 |
| PN-IQ | 333 | 618 | 72 | 1023 |
| PNNN sparse ideal | 179 | 592 | 9 | 780 |

GMP/PN-IQ cuentan 14 raíces por muestra; PN-IQ añade 2 divisiones y una
comparación de cero. PNNN cuenta 14 raíces, 2 divisiones de fase y 12
evaluaciones sigmoides (hasta 12 exponenciales). La aritmética interna de
sigmoide no tiene equivalencia FLOP asignada. No sumar raíces, divisiones
y activaciones como si fueran operaciones homogéneas.

El coste PNNN presupone omitir pesos cero: 296 pesos activos + 14 sesgos.
El coste denso del mismo calendario es 2252. Los campos genéricos
`RealMultiplicationsPerSample`/`CoreFLOPsPerSample` del contador sparse
conservan el baseline denso; para el sparse usar
`IdealSparseRealMultiplicationsPerSample=440`,
`IdealSparseRealAdditionsPerSample=340` y `FLOPsPerSample=780`.
No interpretar los campos densos como el desglose sparse ni sumar dos veces
`NumAbs`/`NumSqrt` o sigmoide/exponencial.

GMP usa el calendario de coeficientes complejos literales y PN-IQ el calendario
real con reutilización estructural existente. Esto no es una comparación
entre implementaciones hardware igualmente optimizadas ni una cota inferior
intrínseca de cada arquitectura. En P310, 620 de los 644 FLOPs ahorrados
corresponden a ponderación/acumulación; se evalúan funciones y soportes distintos.

El calendario no añade explícitamente el z-score de 84 entradas (84 restas,
84 divisiones) ni la transformación afín de las dos salidas (2 productos y
2 sumas). Se pueden fusionar con las capas en otra implementación, pero
no se presupone aquí. Tampoco se cuantifican accesos de memoria, retardos,
DC del registro, transferencias, latencia ni energía. Se preserva el indicador
histórico; no se presenta una cifra ficticia de «coste físico completo».

## Error por amplitud

`amplitude_bins_P310.csv` guarda límites físicos y amplitud normalizada por
el máximo del registro DC-removed (aprox. 0.0665). Son 20 bins por cuantiles
de evaluación, de 22118–22119 muestras cada uno. Se guardan energía de
referencia, energía del error, error medio por muestra, fracción de error
total y NMSE interno del bin. Ningún bin se elimina.

GMP gana a PN-IQ en los siete bins de menor amplitud (35% de las muestras);
PN-IQ gana en los otros trece. La ganancia PN-IQ va desde -0.50 dB en el
primer bin a +0.76 dB en el penúltimo, y +0.4911 dB en el último.
La menor energía de referencia del primer bin (0.118% del total) recomienda
no convertir su NMSE relativo en una afirmación fuerte sobre potencia.

El 5% de muestras de mayor amplitud concentra 32.72% de la energía del
objetivo y 55.07%/54.69%/38.21% del error total de GMP/PN-IQ/PNNN.
El 10% superior concentra 65.99%/64.88%/46.81% del error.
NMSE del último bin: -36.1424/-36.6335/-34.8086 dB, respectivamente.
La menor fracción PNNN en esa cola NO implica menor error absoluto: PNNN
tiene mucho más error global y en amplitudes bajas.

Esto localiza la ventaja PN-IQ en amplitudes medias/altas para esta captura.
No identifica causalmente memoria térmica, load modulation ni otro mecanismo
físico. P100/P500 no se añadieron al análisis de bins para evitar gráficos
redundantes; sus resultados de eficiencia sí están en las tablas completas.

## Figuras y reproducción

Tres figuras inspeccionadas visualmente, cada una en FIG editable, PNG y TikZ:

- [Pareto](../results/dpd_efficiency/20261020_204728/pareto.png): círculos/fronteras; cruces para todos los puntos.
- [Techos FLOPs](../results/dpd_efficiency/20261020_204728/equal_flops.png): escalones de selección discreta, no interpolación experimental.
- [Amplitud P310](../results/dpd_efficiency/20261020_204728/amplitude_error_P310.png): NMSE interno y reparto de energía del error.

`analysis.mat` conserva tablas, predicciones de P310, bins, split,
configuración, procedencia y desglose. La carpeta incluye copia del ejecutor.
Desde la raíz: `addpath('DPD/experiments'); result = analyzeDPDEfficiency;`
crea otra carpeta; para consultar esta ejecución usar `load` sin regenerar.

Conclusión: hay evidencia offline de eficiencia PN-IQ a calidad exigente,
pero GMP conserva ventajas a costes muy bajos y amplitudes bajas. PNNN es
un baseline inferior en esta campaña, no un objetivo de optimización de este
trabajo. Falta confirmar la ventaja con captura independiente y medidas RF.
