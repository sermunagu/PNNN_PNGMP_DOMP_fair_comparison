# DPD offline por imitación del ILC

`testDPD` es ahora una función: `testDPD(P, modo)`, con P=340 y modo
`"reuse"` por defecto. Añade `config/` y `toolbox/` al path. No controla RF.
Usa exclusivamente `measurements/experiment20260429T134032_xy.mat`: `x=u`
deseada, `y=x_ILC` predistorsionada. No usa el MAT `_forward_xy` ni requiere
el repositorio ILC vecino. La firma del MAT coincide con `395e47d7554e…`.

## Estado del experimento existente

Los checkpoints históricos de `results/dpd_ilc/P340_395e47d7554e/` conservan
predicciones, soportes y métricas lineales, pero no coeficientes; `pnnn_P340.mat` tiene máscara y
predicciones, pero no la red sparse final. `pnnn_dense_source.mat` sí guarda
red densa, normalización, firma y duración de fine-tuning. Las predicciones
no permiten reconstruir fielmente los parámetros perdidos.

Recuperación autorizada y completada el 9 de octubre de 2026: se reajustaron
los dos modelos lineales sobre sus soportes guardados (sin DOMP) y se repitió
poda + 164 épocas de fine-tuning desde la red densa existente, sin repetir sus
1228 épocas. La normalización histórica permanece intacta.
Los modelos reutilizables están en `models_P340.mat`, con checkpoints nuevos
`linear_models_P340.mat` y `pnnn_model_P340.mat`. Los cuatro MAT históricos
conservan exactamente su SHA-256; no se sobrescribió ninguno.

## Comandos MATLAB

La recuperación ya ejecutada desde la raíz fue:

```matlab
[bundle, dpd] = testDPD(340, "recover");
```

Para un presupuesto nuevo, identificación explícita de un solo punto (puede
ser costosa; no ejecuta el sweep completo):

```matlab
[bundle, dpd] = testDPD(200, "train");
```

Después de obtener los modelos, la carga no entrena:

```matlab
[bundle, dpd] = testDPD(340); % carga y nueva exportación con timestamp
addpath(genpath('toolbox'));
modelFile = 'results/dpd_ilc/P340_395e47d7554e/models_P340.mat';
dpd = applyDPDModels(modelFile, uNueva, 'mi_nueva_execution.mat');
```

`uNueva` debe ser un registro periódico completo, complejo, finito, más largo
que M, a la misma tasa (491.52 MHz en este experimento), amplitud física y
convención que el entrenamiento. La compatibilidad física de una señal nueva
debe comprobarla el usuario: no se deduce de sus muestras. No es una API de
streaming por bloques; los retardos se extienden circularmente sobre el registro.

La entrada se centra restando su propia media, igual que en identificación.
La salida representa el objetivo ILC sin DC; **no** se restaura su media,
normaliza el pico, recorta ni adapta a límites del DAC. PNNN reutiliza mu/sigma
del entrenamiento y su rotación inversa; PN-IQ conserva su construcción I/Q y
restauración de fase. GMP usa los coeficientes físicos ya ajustados, sin volver
a calcular normalizaciones con la entrada nueva.

## Archivos y controles

Se guardan en la carpeta P/firma del experimento:

- `linear_models_P<P>.mat`: coeficientes físicos y de comparación, definición
  compacta de regresores, soportes, correspondencia I/Q, configuración y métricas.
- `pnnn_model_P<P>.mat`: red **final** sparse, máscaras, normalización, arquitectura,
  métricas y firma de la densa. Ambos sidecars permiten reanudar una interrupción.
- `models_P<P>.mat`: bundle autónomo con los tres modelos, configuración, split,
  medias originales, procedencia, tasa de muestreo, métricas y versión MATLAB.
- `experiment20260429T134032_xy_execution_<timestamp>.mat`: `dpd(k).yvalmod`,
  `dpd(k).modeltype` y metadatos; nunca sustituye el execution original.

Antes de publicar el bundle se comparan predicciones de ajuste e inferencia.
En recuperación se comparan además soportes/máscaras y predicciones antiguas:
error relativo L2 <=1e-10 para lineales (double), <=1e-6 para PNNN (learnables
single). Si no coinciden, se detiene sin aceptar silenciosamente otro modelo.
Los sidecars conservan el trabajo para examinarlo; no borrar los originales.

## Auditoría de P=340

| Familia | Grados reales activos del modelo | Almacenamiento numérico principal |
|---|---|---|
| Complex GMP | 170 coeficientes complejos = 340 reales | 170 complejos físicos; además coeficientes de comparación y metadatos |
| PN-IQ | 170 características, 170 coeficientes I + 170 Q | 340 reales físicos; además comparación y metadatos |
| PNNN | 326 pesos en máscara + 14 biases protegidos | 1046 posiciones densas (1032 pesos + 14 biases), incluidas las podadas a cero |

PNNN usa 84 entradas, fc1(12), sigmoide y fcOut(2); la sigmoide no tiene
parámetros entrenables. Activo significa grado permitido por soporte/máscara,
no `nnz` del valor numérico: un bias activo puede valer cero. No se cuentan
estadísticas mu/sigma, regresores, índices, máscaras ni constantes de activación
como parámetros ajustados. El criterio compara grados ajustables, no bytes de
almacenamiento ni capacidad funcional idéntica. Las columnas redundantes o
estructuralmente constantes pueden reducir grados efectivos: no se altera aquí
el criterio histórico ni se elimina ninguna característica.

FLOPs es otra magnitud: el checkpoint reporta 1817/1107/840 por muestra,
respectivamente. La cifra sparse de PNNN supone aprovechar pesos podados;
MATLAB almacena matrices densas y no garantiza ese coste real de ejecución.

Limitación numérica detectada en la prueba sintética: el canal imaginario del
tap actual, teóricamente cero tras normalizar fase, puede tener sigmaX~1e-17.
La regla histórica solo sustituye sigma exactamente cero; puede amplificar
redondeo y romper la equivariancia de fase numérica (1.03% L2 en este fixture).
La inferencia conserva esa regla y reproduce la evaluación original. Una
corrección futura sería fijar explícitamente ese canal a cero o umbralizar su
sigma, pero requiere aprobación y nueva validación/entrenamiento; no se aplica
silenciosamente ni se atribuye ese porcentaje al experimento real.
En el checkpoint denso real se confirma sigmaX=8.019388e-19 en la característica
15; no se ha medido aquí el efecto de rotar la captura real.

## Validación y límites

```matlab
addpath('tests'); run_dpd_model_reuse_test
run('tests/run_linear_complexity_sweep_test.m')
run('tests/run_pnnn_shared_dense_sweep_test.m')
```

El test nuevo es sintético (384 muestras, una época sparse, sin entrenamiento
denso): comprueba P340, serialización, nuevas entradas, memoria periódica,
rotación, formato de exportación, protección de archivos y recuperación de
soportes. No confundir ese test con validación experimental RF.

La validación real posterior a la recuperación también pasó en otra sesión de
MATLAB: error relativo L2 exactamente 0 para los tres modelos frente al ajuste
recuperado, las predicciones históricas y la nueva exportación. Tolerancias:
1e-10 para lineales double y 1e-6 para PNNN con learnables single. Coinciden
soportes, máscaras, normalización, presupuesto, lambdas y FLOPs. Informe:
`validation_20261009_001045_626.mat`; señales nuevas:
`experiment20260429T134032_xy_execution_20261009_000954_626.mat`.

Validación ejecutada en MATLAB R2026a: las tres pruebas anteriores pasan;
`checkcode` no emite diagnósticos en los nueve archivos MATLAB modificados o
creados. El modo por defecto se detiene correctamente si faltan modelos,
sin entrenar ni exportar. Los cuatro MAT históricos conservan su SHA-256.

NMSE histórico offline: Complex −38.456 dB, PN-IQ −38.931 dB, PNNN −35.630 dB,
respecto a x_ILC centrada, no a la salida del PA. Quedan pendientes transmisión
autorizada, límites/amplitud del DAC, decisión sobre DC, sincronización, y medida
de NMSE/ACLR/espectro de salida del PA; este flujo no realiza ninguna de ellas.
