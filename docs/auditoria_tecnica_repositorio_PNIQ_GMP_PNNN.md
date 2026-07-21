# Auditoría técnica exhaustiva del pipeline GMP–PN-IQ-GMP–PNNN

**Repositorio:** `sermunagu/PNNN_PNGMP_DOMP_fair_comparison`  
**Rama auditada:** `main`  
**Commit auditado:** `675687dc9b5a4b5c603a533f8b39ef5a53413eac`  
**Finalidad:** documento interno previo a la redacción del primer paper. No es todavía el texto del congreso.

## 1. Alcance real de la auditoría

Esta auditoría cubre de forma exhaustiva el **pipeline científico canónico ejecutado por `main_sweep_and_comparison`**, sus modelos, selección de soporte, estimación, pruning, métricas, complejidad, selección del punto operativo, figuras, checkpoints y pruebas directamente relacionadas. También cubre el experimento auxiliar widely-linear.

No se ha auditado línea por línea el código externo de `third_party/matlab2tikz`, porque es un submódulo de terceros. Tampoco puede inspeccionarse semánticamente el contenido binario completo de los `.mat` mediante el conector de GitHub; sí se auditan todos los contratos mediante los que el código los carga, valida y utiliza.

## 2. Conclusión ejecutiva

El repositorio implementa tres familias con presupuestos **exactos de parámetros reales**:

1. **Complex GMP-DOMP:** para un presupuesto \(P\), selecciona \(P/2\) regresores complejos y ajusta \(P/2\) coeficientes complejos.
2. **PN-IQ-GMP con selección DOMP:** para \(P\), selecciona \(P/2\) características reales de un diccionario I/Q normalizado en fase y ajusta dos vectores reales de \(P/2\) coeficientes, uno por salida I y otro por salida Q.
3. **Sparse PNNN N12:** poda una red de una capa oculta con 12 neuronas hasta exactamente \(P\) escalares activos, protegiendo siempre los 14 sesgos.

La implementación principal es científicamente coherente y contiene controles fuertes de reproducibilidad. Los puntos que deben cerrarse antes del paper son:

- quitar `clear` del interior de `main_sweep_and_comparison`, porque elimina variables del workspace de la función y puede romper el argumento manual;
- abandonar el nombre **PN-DOMP**: DOMP no está normalizado en fase; DOMP se aplica al diccionario PN-IQ;
- no llamar “test independiente” a la evaluación sobre la señal completa, porque incluye las muestras de identificación;
- dejar explícito que los FLOPs de la PNNN dispersa suponen un kernel que omite pesos cero;
- no interpretar cuantitativamente la curva de máximo coeficiente de PNNN frente a las lineales, porque usan parametrizaciones distintas;
- explicar que igualdad de parámetros no implica igualdad de espacio candidato.

## 3. Mapa del flujo canónico

`main_sweep_and_comparison.m`:

1. carga la configuración;
2. ejecuta `run_parameter_sweep`;
3. obtiene el primer presupuesto conjuntamente estabilizado;
4. permite un override manual;
5. ejecuta `run_selected_comparison`.

`run_parameter_sweep.m`:

1. verifica la herramienta de exportación;
2. carga `x,y`;
3. aplica el mapeo local `xy_forward`;
4. elimina DC;
5. crea la identificación común y el split interno;
6. construye o reutiliza el barrido lineal;
7. construye o reutiliza las referencias Ridge;
8. construye o reutiliza una única PNNN densa N12;
9. poda independientemente esa fuente para cada presupuesto;
10. une las tres tablas;
11. selecciona el punto operativo;
12. exporta CSV, MAT y figuras.

## 4. Datos, mapeo y particiones

La configuración usa:

- `mappingMode = 'xy_forward'`;
- fracción de identificación 0.10;
- semilla de identificación 1004;
- 85 % de la identificación para entrenamiento interno;
- semilla interna 42;
- 10 bins de amplitud;
- eliminación de DC activa.

La identificación **no es un muestreo aleatorio del 10 %**. `sel_indices` toma un bloque contiguo de aproximadamente el 10 % que contiene el máximo de \(|y|\). El split train/validation equilibrado en amplitud se conserva exclusivamente para seleccionar el entrenamiento y fine-tuning de la PNNN; los modelos lineales usan toda la identificación.

La señal completa son todas las muestras. Por tanto, incluye la identificación y no constituye un holdout independiente.

## 5. Población GMP

El manager inicializa una población determinista GMP con:

- orden máximo 13;
- memoria base `L=10`;
- offset cruzado `M=2`;
- límites generales `Qpmax=Qnmax=50`;
- evaluación genética desactivada para el pipeline actual: se usa la población inicial completa y después DOMP.

Los términos canónicos implementados son:

\[
u^{(a)}_{k,\ell}[n]=x[n-\ell]|x[n-\ell]|^k,
\quad k=0,\ldots,12,\quad \ell=0,\ldots,10,
\]

\[
u^{(b)}_{k,\ell,m}[n]=x[n-\ell]|x[n-\ell-m]|^k,
\quad k=1,\ldots,12,\quad m=1,2,
\]

\[
u^{(c)}_{k,\ell,m}[n]=x[n-\ell]|x[n-\ell+m]|^k.
\]

Se incluyen todos los grados totales 1–13, no únicamente órdenes impares.

Además, `removerepeated` fuerza inicialmente tres auxiliares: \(x[n]\), \(x^*[n]\) y \(|x[n]|\); después elimina duplicados. El término \(x[n]\) ya aparece en el GMP y se deduplica, mientras los otros dos quedan como auxiliares reconocidos. El tamaño resultante esperado es 673 columnas complejas: 671 canónicas y dos auxiliares.

La evaluación de regresores usa indexado periódico:

\[
x[n-q]\equiv x[((n-q-1)\bmod N)+1].
\]

## 6. DOMP exacto

`selectDOMPSupport(X,y,K,\tau)` no ajusta Ridge durante la selección.

Inicialización:

\[
Z^{(0)}=X,\qquad r^{(0)}=y,\qquad \mathcal S_0=\varnothing.
\]

En la iteración \(t\):

1. se calculan las normas de las columnas residuales de \(Z\);
2. se descartan columnas ya seleccionadas o por debajo de la tolerancia;
3. las columnas elegibles se normalizan;
4. se calcula
   \[
   s_j=\left\|z_j^H r^{(t-1)}\right\|_2;
   \]
5. se selecciona el máximo;
6. se ortogonalizan todos los candidatos restantes:
   \[
   Z\leftarrow Z-q(q^H Z);
   \]
7. se recalculan coeficientes y residuo sobre la matriz original:
   \[
   \hat c_{\mathcal S_t}
   =\arg\min_c\|y-X_{\mathcal S_t}c\|_2^2,
   \quad
   r^{(t)}=y-X_{\mathcal S_t}\hat c_{\mathcal S_t}.
   \]

El LS usa QR si la matriz se considera de rango completo y `lsqminnorm` si detecta rango deficiente.

Cuando \(X\) es real y \(y\) complejo —caso PN-IQ— la correlación es compleja y su módulo combina simultáneamente la información de I y Q.

## 7. Complex GMP-DOMP

### 7.1 Ruta de identificación

Se construye una matriz compleja GMP sobre toda la identificación. Se calcula **una única ruta DOMP máxima de 250 regresores complejos**, correspondiente a 500 parámetros reales. Cada presupuesto usa un prefijo anidado:

\[
K=P/2.
\]

### 7.2 Ajuste principal sin regularización

Para cada prefijo, se normalizan columnas:

\[
D=\operatorname{diag}(\|u_1\|_2,\ldots,\|u_K\|_2),\qquad
\widetilde U=UD^{-1}.
\]

Para \(\lambda=0\):

\[
\widetilde h=\widetilde U^\dagger y
\]

mediante `lsqminnorm` y una tolerancia de rango robusta.

Los coeficientes físicos son:

\[
h=D^{-1}\widetilde h.
\]

La solución principal usa siempre \(\lambda=0\). `cfg.lambdaGrid` se conserva solo para experimentos auxiliares y no participa en el pipeline lineal canónico.

### 7.3 Evaluación

Cada prefijo se ajusta una sola vez sobre identificación y se evalúa sobre identificación y señal completa. La señal completa se procesa por bloques de 8192 muestras.

### 7.4 Parámetros

Con \(K\) coeficientes complejos:

\[
P_{\mathrm{GMP}}=2K.
\]

## 8. PN-IQ-GMP propuesto

### 8.1 Rotación

Para cada muestra:

\[
r[n]=
\begin{cases}
x^*[n]/|x[n]|,&|x[n]|\neq0,\\
1,&|x[n]|=0.
\end{cases}
\]

Se rota el objetivo:

\[
z[n]=r[n]y[n].
\]

Cada regresor GMP se rota:

\[
\widetilde u_j[n]=r[n]u_j[n].
\]

### 8.2 Diccionario real

De cada regresor se crean:

\[
f_{j,I}[n]=\Re\{\widetilde u_j[n]\},\qquad
f_{j,Q}[n]=\Im\{\widetilde u_j[n]\}.
\]

Para regresores canónicos cuyo portador es \(x[n]\), la parte Q es estructuralmente cero y se elimina. Los auxiliares se calculan mediante la rotación exacta.

Partiendo de 673 regresores complejos, el diccionario PN-IQ esperado tiene 1285 columnas reales, porque 61 componentes Q canónicas con portador de retardo cero son nulas.

### 8.3 DOMP

DOMP se aplica a:

- matriz candidata **real** \(F\);
- objetivo complejo rotado \(z\).

No existe un algoritmo distinto llamado PN-DOMP. Es DOMP estándar aplicado al diccionario PN-IQ.

Se construye una ruta máxima de 250 características reales. Para un presupuesto \(P\):

\[
M=P/2
\]

características reales activas.

### 8.4 Ajuste I/Q independiente

Sobre el mismo soporte \(\mathcal S\):

\[
\hat z_I=F_{\mathcal S}c_I,
\qquad
\hat z_Q=F_{\mathcal S}c_Q.
\]

Las dos salidas comparten soporte, se ajustan independientemente mediante `lsqminnorm` y usan siempre \(\lambda=0\) en la solución principal:

\[
\hat z=F_{\mathcal S}c_I+jF_{\mathcal S}c_Q,
\qquad
\hat y[n]=r^*[n]\hat z[n].
\]

### 8.5 Libertad adicional

Para \(\widetilde u=a+jb\), un coeficiente complejo convencional \(h\) impone:

\[
h\widetilde u=h\,a+(jh)\,b.
\]

Por tanto, los coeficientes complejos asociados a \(a\) y \(b\) están restringidos por:

\[
g_b=jg_a.
\]

PN-IQ permite asignar coeficientes complejos independientes a cada característica real seleccionada:

\[
\hat z=g_a a+g_b b,
\]

sin imponer \(g_b=jg_a\). Además, DOMP puede seleccionar I y Q de un mismo regresor por separado.

La formulación es más flexible, pero debe explicarse que su espacio candidato tiene más átomos escalares que el GMP complejo. Igual número de parámetros activos no significa igual número de candidatos.

### 8.6 Parámetros

Con \(M\) características reales y dos coeficientes por característica:

\[
P_{\mathrm{PN-IQ}}=2M.
\]

## 9. Referencias Fixed Ridge

Las variantes Ridge:

- no recalculan DOMP;
- reutilizan exactamente los soportes principales;
- prueban \(\lambda=10^{-3},10^{-4},10^{-5}\);
- conservan parámetros y FLOPs;
- cambian únicamente la estimación de coeficientes.

Sirven para mostrar el compromiso entre NMSE y rango dinámico numérico.

## 10. PNNN

### 10.1 Vector de entrada

Con \(M=13\) y órdenes \(\{1,3,5,7\}\), para cada muestra se construyen 14 taps periódicos:

\[
x_n=[x[n],x[n-1],\ldots,x[n-13]].
\]

Se aplica la misma rotación \(r[n]\). En modo `full`:

\[
\phi[n]=[
\Re\{r[n]x_n\},
\Im\{r[n]x_n\},
|x_n|^1,
|x_n|^3,
|x_n|^5,
|x_n|^7]^T.
\]

La dimensión es:

\[
D=2(14)+4(14)=84.
\]

El objetivo es:

\[
t[n]=[
\Re\{r[n]y[n]\},
\Im\{r[n]y[n]\}]^T.
\]

Observación de auditoría: en modo full, \(\Im\{r[n]x[n]\}=0\) estructuralmente, y \(\Re\{r[n]x[n]\}=|x[n]|\), que duplica la característica de envolvente de orden 1 para el tap actual.

### 10.2 Normalización

Se usa z-score por característica y por canal de salida:

\[
\bar\phi_d=(\phi_d-\mu_{X,d})/\sigma_{X,d},
\qquad
\bar t_q=(t_q-\mu_{Y,q})/\sigma_{Y,q}.
\]

Las desviaciones cero se sustituyen por 1.

Durante selección, las estadísticas se calculan solo en train interno. En el refit final se recalculan sobre toda identificación.

### 10.3 Arquitectura

\[
h=\sigma(W_1\bar\phi+b_1),
\qquad
\bar{\hat t}=W_2h+b_2,
\]

con:

- entrada 84;
- 12 neuronas ocultas;
- activación sigmoide;
- salida lineal de 2 canales.

Parámetros densos:

\[
P_{\mathrm{dense}}
=84\cdot12+12+12\cdot2+2
=1046.
\]

### 10.4 Entrenamiento denso

Se usa Adam, MSE, batch 1024, learning rate inicial \(2\cdot10^{-4}\), schedule piecewise, factor 0.95, shuffle por época y CPU.

El presupuesto de actualizaciones se escala para igualar un protocolo histórico basado en 70 % de la señal y 150 épocas. La red de selección se entrena en train interno y se conserva la de menor validation loss. Después se reinicializa con la misma semilla y se vuelve a entrenar desde cero sobre toda identificación durante el número de épocas seleccionado.

### 10.5 Selección del número de épocas de fine-tuning

Se poda la red de selección al presupuesto de referencia 200. Sobre ese modelo se ejecuta fine-tuning con validation y se elige una única duración. Esa duración se reutiliza para todos los presupuestos finales.

### 10.6 Poda

La poda es global por magnitud sobre todos los pesos de `fc1` y `fcOut`.

- Los sesgos nunca se podan.
- Hay 14 sesgos protegidos.
- Para un objetivo \(P\), quedan \(P-14\) pesos y 14 sesgos.
- Las máscaras se derivan siempre de la misma fuente densa inmutable.
- Cada presupuesto se fine-tunea independientemente.
- Los gradientes podados se enmascaran.
- Después de cada actualización Adam se reaplica la máscara.
- Se verifica que ningún peso podado vuelva a ser no nulo.

### 10.7 Predicción

La red produce dos canales normalizados. Se desnormalizan, se forma:

\[
\hat z=\hat t_I+j\hat t_Q,
\]

y se restaura:

\[
\hat y[n]=r^*[n]\hat z[n].
\]

## 11. NMSE

Todos los modelos usan:

\[
\operatorname{NMSE}_{\mathrm{dB}}
=10\log_{10}
\frac{\sum_n|y[n]-\hat y[n]|^2}
{\sum_n|y[n]|^2}.
\]

La misma función se usa para identificación y señal completa.

## 12. FLOPs

Convención:

- suma real: 1;
- multiplicación real: 1;
- suma compleja: 2;
- multiplicación compleja: 6;
- MAC complejo: 8.

Raíces, divisiones, valores absolutos y activaciones se cuentan aparte, pero no se convierten a FLOPs. Por tanto, `FLOPsPerSample` es complejidad aritmética básica, no coste total de hardware ni tiempo de ejecución.

Para PNNN dispersa, los FLOPs reportados asumen un kernel ideal que omite pesos exactamente cero. MATLAB sigue evaluando matrices densas durante `predict`; el ahorro es una estimación de implementación sparse.

## 13. Máximo coeficiente

### 13.1 Modelos lineales

La definición solicitada por el tutor se implementa así:

1. entrada y salida con máximo módulo uno en identificación;
2. construcción del diccionario;
3. normalización L2 de cada columna;
4. ajuste en esa base;
5. máximo absoluto de los escalares reales activos.

Para GMP complejo:

\[
C_{\max}=
\max_k\{
|\Re\{\widetilde h_k\}|,
|\Im\{\widetilde h_k\}|
\}.
\]

Para PN-IQ:

\[
C_{\max}=
\max_k\{
|(\widetilde c_I)_k|,
|(\widetilde c_Q)_k|
\}.
\]

El código usa las señales originales y divide los coeficientes de columnas normalizadas por el pico de salida. La normalización global de entrada cancela al normalizar cada columna. La prueba `run_linear_complexity_sweep_test` reconstruye explícitamente la normalización de pico y verifica igualdad numérica.

### 13.2 PNNN

La PNNN usa el máximo de todos los pesos y sesgos activos en su parametrización z-score. No es la misma base que las columnas GMP normalizadas. La figura es válida como diagnóstico de rango numérico por implementación, pero no como comparación directa de ganancias físicas.

## 14. Punto operativo

El punto se selecciona buscando el menor presupuesto común para el que ninguna de las tres familias mejora más de 0.20 dB dentro de los 100 parámetros siguientes.

El primer punto es 340 parámetros:

- ganancia futura GMP: 0.0112 dB;
- ganancia futura PN-IQ: 0.1952 dB;
- ganancia futura PNNN: 0.0772 dB.

## 15. Resultados en 340 parámetros

| Familia | NMSE señal completa | FLOPs/muestra | Máx. parámetro |
|---|---:|---:|---:|
| Complex GMP-DOMP | -38.85998 dB | 1807 | 6.34896e6 |
| PN-IQ-GMP, DOMP | -42.74017 dB | 1095 | 4.59256e6 |
| Sparse PNNN N12 | -39.24296 dB | 840 | 2.17243 |

PN-IQ mejora 3.88019 dB frente a GMP y reduce FLOPs.

Las variantes Ridge en 340 reducen el máximo coeficiente a aproximadamente 86–132, a cambio de empeorar NMSE.

## 16. Espectro y tiempo

Welch común:

- ventana Hann periódica 2048;
- solapamiento 75 %;
- NFFT 16384;
- bilateral centrado;
- error \(e=y-\hat y\);
- referencia: máximo de la PSD de la salida medida.

Las figuras temporales usan las 1000 muestras centrales, con magnitud y parte real o imaginaria. La superposición visual no contradice el NMSE: las diferencias viven en el residuo y son pequeñas frente a la señal principal.

## 17. Probe widely-linear

`run_widely_linear_gmp_probe` concatena:

\[
[U,U^*]
\]

y selecciona 170 átomos para mantener 340 parámetros reales. Resultado:

- GMP: -38.85998 dB;
- WL-GMP: -38.86672 dB;
- PN-IQ: -42.74017 dB.

La mejora WL es solo 0.00674 dB. Esto indica que añadir \(U^*\) con coeficientes fijos no reproduce la estructura PN-IQ. El probe es auxiliar y no entra en el sweep principal.

## 18. Checkpoints y reproducibilidad

La identidad del sweep incluye:

- firma de la señal y configuración;
- rejilla de parámetros, lambdas Ridge fijas y protocolo lineal;
- población GMP;
- opciones DOMP;
- configuración PNNN;
- configuración de entrenamiento y pruning.

Se usa SHA-256. Cada checkpoint se rechaza si pertenece a otra identidad. La fuente densa PNNN tiene además una firma de learnables y estadísticas de normalización.

## 19. Hallazgos de auditoría

### Crítico / corregir antes de usar override manual

`main_sweep_and_comparison` contiene:

```matlab
clear; clc; close all force;
```

`clear` elimina todas las variables del workspace actual, incluidos los argumentos de la función. La ejecución sin argumento funciona por cortocircuito de `nargin`, pero el override manual puede quedar sin `selectedParameters`. Debe eliminarse `clear`.

### Alto

1. **Nombre PN-DOMP:** incorrecto conceptualmente. Usar `PN-IQ-GMP` y “DOMP-based support selection”.
2. **Full signal no independiente:** nunca denominarlo test independiente.
3. **FLOPs sparse ideales:** declararlo explícitamente.
4. **Máximo coeficiente PNNN:** parametrización distinta.
5. **Espacios candidatos distintos:** la comparación es por parámetros activos, no por tamaño de diccionario.

### Medio

1. El modo PNNN `full` contiene una característica cero y una duplicada.
2. `run_widely_linear_gmp_probe` fija directamente el directorio `sweep_d113e389ab78`; es frágil si cambia la identidad.
3. El nombre de la tabla y las leyendas siguen usando `PN-IQ PN-DOMP`.
4. El título “Max. abs. real coefficient” sería más preciso como “maximum absolute active real scalar component”.
5. `exportPaperFigure` restaura `Visible` después de `savefig`, pero no mediante `onCleanup`; un error durante `savefig` puede dejar la figura visible.

## 20. Respuestas que debes dominar

- Cada familia lineal ejecuta DOMP una sola vez sobre toda identificación.
- Cada presupuesto usa un prefijo de esa ruta y se ajusta con \(\lambda=0\).
- DOMP no se aplica a coeficientes: se aplica a columnas candidatas del diccionario.
- GMP usa un diccionario complejo y objetivo complejo.
- PN-IQ usa un diccionario real I/Q y objetivo complejo rotado.
- En PN-IQ el soporte es común para las dos salidas; los coeficientes I y Q se estiman de forma independiente.
- Ridge no cambia el soporte.
- La PNNN no usa DOMP.
- La PNNN se poda globalmente por magnitud, protege sesgos y parte siempre de la misma red densa.
- Los 340 parámetros incluyen todos los pesos y sesgos activos.
- El full-signal NMSE incluye identificación.
- Los FLOPs PNNN son ideales bajo un kernel sparse.
- El máximo coeficiente lineal sí sigue la normalización solicitada; la PNNN no está en la misma parametrización.
