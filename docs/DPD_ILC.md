# DPD offline por imitación del ILC

## OFFLINE / ONLINE

Desde la raíz del repositorio:

```matlab
modelFile = trainDPDOffline('measurements/nueva_medida.mat', 340);
outputFile = runDPDOnline(modelFile, 'measurements/entrada_online.mat');
```

OFFLINE carga directamente `x` e `y`. Cada llamada ajusta de nuevo los tres
modelos, conservando el centrado DC, la identificación, DOMP, N12/sigmoide,
selección de épocas, poda y fine-tuning existentes. No reutiliza entrenamientos.

Una carpeta nueva creada con `tempname` bajo `results/dpd_offline/` contiene:

- `models.mat`: los tres modelos en `bundle`, configuración, split, métricas,
  archivo de origen y `fs` si existe.
- `predictions.mat`: predicciones directas del ajuste.

ONLINE carga el `bundle`, lee `x` y llama a `applyDPDModels`, sin entrenar.
Guarda un MAT nuevo junto al modelo, con `dpd(k).yvalmod` y `dpd(k).modeltype`
en orden GMP-DOMP, PN-IQ-GMP, PNNN. Tras reiniciar MATLAB, usar la ruta del
`models.mat` guardado. Los paquetes existentes siguen siendo compatibles.

Las rutas relativas las resuelve MATLAB. Se usa el formato conocido del tutor;
archivo o campos incorrectos producen los errores normales de MATLAB.
No hay helpers de carga, hashes de medidas, UUID ni comprobaciones de calibración.
El investigador debe usar señales compatibles en frecuencia, escala y periodicidad.
No se escala, recorta ni restaura DC; el NMSE offline no es una medida RF.

Se conservan los tests sintéticos ligeros de entrenamiento y recarga. El tercer
argumento opcional `cfg` de OFFLINE permite reducir épocas y usar una carpeta
temporal en esos tests; las dos llamadas normales mantienen el protocolo original.

Tras la simplificación pasan los tests sintéticos y la recarga en otra sesión
MATLAB. El P340 existente reproduce exactamente, señal por señal y con el mismo
tipo numérico, su archivo ONLINE anterior (error relativo 0 en los tres modelos).
No se repitió su entrenamiento. `checkcode` no informa diagnósticos.

## Barrido de modelado PA y reanudación

El barrido canónico sigue usando por defecto la captura `_forward_xy.mat`;
no es el futuro barrido DPD ni cambia el significado físico de x/y.

```matlab
sweep = run_parameter_sweep(); % nueva carpeta, no reutiliza otros experimentos
sweep = run_parameter_sweep('results/parameter_sweep/carpeta_de_esta_ejecucion');
comparison = run_selected_comparison(sweep.selection.selectedParameters, sweep);
```

La segunda llamada es una alternativa para reanudar, no otro entrenamiento.
La carpeta se imprime antes de entrenar. `experiment.mat` guarda configuración,
señales preprocesadas y split; los MAT lineal, Ridge, denso y
`pnnn_target_XXXX.mat` guardan cada etapa. Solo se calculan etapas ausentes.
Cada PNNN parte de la misma fuente densa, no de la poda anterior.
No cambie los archivos de esa carpeta para reutilizarla con otra captura.

Los nuevos checkpoints usan variables directas (`linear`, `fixedResults`,
`denseSource`, `point`, `result`), sin firmas ni versiones de protocolo.
No se migran ni sobrescriben checkpoints históricos: para ejecutarlos con el
formato anterior debe usarse la revisión histórica del código.
Para regenerar figuras de una ejecución nueva se carga
`complexity_sweep.mat` y se pasa `S.result` a `run_selected_comparison`.
Reanudar vuelve a exportar las tres figuras del barrido, sin repetir ajustes.

Las métricas y checkpoints se guardan antes de exportar. No hay preflight:
`matlab2tikz` solo se necesita al producir figuras FIG/PNG/TikZ. Un error gráfico
no borra el entrenamiento; se puede reanudar la misma carpeta.
El tercer argumento opcional `cfg` permite pruebas sintéticas reducidas,
como en OFFLINE; el uso normal mantiene el preset científico.

## Experimentos históricos

La recuperación P340 del 9 de octubre reajustó los modelos lineales sobre
soportes guardados y repitió únicamente poda/fine-tuning desde la densa.
No se repitieron sus 1228 épocas. Los archivos históricos no se migran.
`testDPD.m` ya estaba eliminado del working tree al iniciar esta limpieza;
no es una llamada operativa vigente. Para inferir use el `models.mat` existente
con ONLINE. La captura DPD es `experiment20260429T134032_xy.mat`:
`x=u`, `y=x_ILC`; no debe confundirse con el modelado PA.

El registro de entrada es periódico. ONLINE centra su propia media y conserva
la normalización entrenada. No restaura DC ni adapta pico, CFR o límites DAC.
La compatibilidad física de las nuevas señales debe comprobarla el investigador.

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
addpath('tests');
run_dpd_offline_online_test;
run_dpd_model_reuse_test;
run_pnnn_shared_dense_sweep_test;
run_sweep_resume_test;
```

Los tests usan capturas sintéticas y pocas épocas. Comprueban soportes,
presupuestos, máscaras, guardado/recarga, ONLINE sin entrenamiento y
reanudación completa o con un único punto ausente. El test de regresión de
simplificación permite registrar una referencia antes de editar y compararla
después; no genera su propio resultado esperado durante la comprobación.

Siguen pendientes la evaluación sobre muestras externas, la decisión sobre
sigma casi nula, y la transmisión y medida RF autorizadas. El NMSE offline
sobre señal completa incluye muestras de identificación; no es un holdout
independiente ni una medida de linealización del PA.
