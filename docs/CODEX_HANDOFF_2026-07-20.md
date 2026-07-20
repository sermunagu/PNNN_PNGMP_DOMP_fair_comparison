# Traspaso de contexto: barrido, figuras de artículo y selección del punto operativo

Fecha del snapshot: 2026-07-20 (Europe/Madrid)

Este documento permite continuar el trabajo en otro chat sin depender del historial de la conversación anterior. Describe por separado:

1. el estado versionado y reproducible en `HEAD`;
2. la implementación científica y técnica ya realizada;
3. las validaciones conocidas;
4. el trabajo no confirmado que estaba siendo modificado simultáneamente por otro chat/proceso.

> **Regla crítica de continuación:** antes de tocar nada, volver a ejecutar las comprobaciones Git de la sección siguiente. El repositorio estaba siendo modificado simultáneamente mientras se escribió este documento. No restaurar, mover, borrar, incluir en un commit ni sobrescribir cambios sin identificar primero su propietario y propósito.

## 1. Repositorio y estado de referencia

- Repositorio local: `C:\Sergi\Investigacion\Códigos\NN\PNNN_PNGMP_DOMP_fair_comparison`
- Remoto esperado: `https://github.com/sermunagu/PNNN_PNGMP_DOMP_fair_comparison.git`
- Rama observada: `feature/paper-figures-selection-and-coefficient-range`
- Seguimiento observado: `origin/feature/paper-figures-selection-and-coefficient-range`
- `HEAD` observado al preparar este traspaso: `c2e4bf198a94deccc9b4161c3d96dac398f7232c`
- Submódulo `third_party/matlab2tikz` observado en: `806c97d99f87f8a1e99a7c54e853c25c82aac301`
- Rama y remoto estaban sincronizados en ese instante (`ahead 0`, `behind 0`).

Comprobación obligatoria al abrir un chat nuevo:

```powershell
Set-Location 'C:\Sergi\Investigacion\Códigos\NN\PNNN_PNGMP_DOMP_fair_comparison'
Get-Location
git rev-parse --show-toplevel
git branch --show-current
git rev-parse HEAD
git status --short --branch
git remote -v
git submodule status
git log --oneline --decorate -8
```

No asumir que `HEAD` o el working tree siguen iguales después de la fecha de este documento.

## 2. Historial útil y fuentes de verdad

Commits relevantes observados, de más reciente a más antiguo:

```text
c2e4bf1 Store validated schema-v3 parameter sweep
840dd51 Add paper figures and quantified operating-point selection
4b649b4 Add paper figures, operating-point selection, and coefficient-range analysis
29492e7 Store validated sweep from current code
```

El commit `29492e7` corresponde al punto de partida estable anterior a la ampliación de figuras, selección automática y métrica de rango. Los commits `4b649b4` y `840dd51` contienen la implementación. `c2e4bf1` incorpora el barrido completo validado con esquema 3 y sus artefactos.

Orden recomendado de autoridad para resolver discrepancias:

1. código y metadatos del `HEAD` actual;
2. artefactos versionados del barrido `sweep_d113e389ab78`;
3. tests del repositorio;
4. este documento como explicación de intención;
5. cambios no confirmados del working tree, solo después de auditarlos.

## 3. Objetivo científico original

La ampliación debía producir material reproducible para artículo y automatizar una selección defendible del presupuesto de PNNN, manteniendo intacta la comparación científica existente. Los objetivos fueron:

- añadir a los barridos el máximo valor absoluto de los parámetros activos de cada modelo;
- seleccionar automáticamente un punto de operación PNNN casi óptimo y de complejidad mínima;
- generar figuras publicables con una paleta IEEE consistente;
- exportar cada figura desde el mismo handle MATLAB a `.fig`, `.png`, `.pdf` y `.tikz`;
- generar espectros solamente para el punto seleccionado;
- validar de antemano que `matlab2tikz` y LaTeX funcionan;
- conservar el split, los objetivos, la definición DOMP, el cálculo de NMSE, el entrenamiento y el significado científico de los modelos.

No se debía cambiar silenciosamente ninguna decisión científica para obtener una figura más favorable. Los cambios autorizados se limitaban a nueva instrumentación, selección cuantificada, presentación y exportación.

## 4. Implementación realizada

### 4.1 Esquema de barrido y firma de experimento

La configuración de barrido se elevó a `schemaVersion = 3`. La versión del esquema forma parte de la firma del experimento para impedir que resultados incompatibles se reutilicen como si fueran equivalentes.

Las preferencias puramente de presentación y selección no se incorporaron a la identidad científica del experimento cuando no afectaban al entrenamiento ni a los datos. Esto permite regenerar presentación o reevaluar un criterio sin fingir que se ejecutó un experimento científico distinto.

Archivos principales:

- `config/getFairDOMPComparisonConfig.m`
- `toolbox/sweep/buildExperimentSignature.m`
- `run_parameter_sweep.m`

### 4.2 Métrica `MaxAbsRealParameter`

Se añadió `MaxAbsRealParameter` a las tablas del barrido principal y del barrido de lambdas fijas.

Definición exacta de intención:

> Maximum absolute active real scalar in the model's native stored parameterization.

Interpretación y cautela:

> The metric compares stored numerical dynamic range for implementation and quantization. The parameters belong to different model parameterizations, so their magnitudes must not be interpreted as identical physical gains.

La métrica se calcula en la parametrización almacenada de cada familia:

- **GMP complejo:** máximo entre los valores absolutos de las partes real e imaginaria de los coeficientes activos.
- **PN:** máximo absoluto de los coeficientes I/Q almacenados.
- **PNNN dispersa:** máximo absoluto de los parámetros aprendibles finales tras el fine-tuning, contando únicamente pesos y sesgos activos según la máscara.
- **Ridge fijo:** mismo concepto aplicado a los coeficientes del modelo correspondiente.

La métrica sirve para discutir rango numérico, representación o cuantización. No convierte parámetros de arquitecturas distintas en ganancias físicas directamente comparables.

Las tablas de esquema 3 quedaron con 13 columnas en el barrido principal y 8 columnas en el barrido de lambda fija. Los lectores/cargadores validan el esquema para evitar mezclar tablas anteriores.

### 4.3 Selección automática del punto operativo

La lógica está centralizada en:

- `toolbox/sweep/selectOperatingPoint.m`

El criterio reportado es:

```text
near-optimal minimum-complexity criterion
```

Para cada presupuesto PNNN:

1. se toma la mejor NMSE PNNN obtenida en ese presupuesto;
2. se compara con GMP al mismo presupuesto;
3. se calcula la pérdida respecto de la mejor NMSE PNNN global;
4. se marca como casi óptimo si la pérdida no excede la tolerancia;
5. se exige además que PNNN no sea peor que GMP al mismo presupuesto;
6. entre los puntos admisibles se escoge el de menor FLOP;
7. si hay empate, se escoge el de menor número real de parámetros.

Tolerancia predeterminada: `0.20 dB`.

Sensibilidades configuradas: `0.10`, `0.15`, `0.20` y `0.25 dB`.

La selección produce diagnósticos completos, no solo el índice elegido: presupuesto, NMSE de PNNN, mejor NMSE global, pérdida frente al mejor, NMSE de GMP emparejado, ganancia frente a GMP, FLOP, parámetros reales, condición de casi optimalidad, condición de superioridad frente a GMP y admisibilidad final.

Artefactos:

- `operating_point_selection.csv`
- `operating_point_selection_sensitivity.csv`
- `operating_point_selection_summary.txt`

La regresión histórica de esquema 2 y el barrido validado de esquema 3 seleccionan el presupuesto `340`. En la regresión histórica conocida:

- presupuesto seleccionado: `340`;
- FLOP seleccionados: `840`;
- NMSE PNNN seleccionada: `-39.242963543132 dB`;
- mejor NMSE PNNN global: `-39.434236892018 dB`, en presupuesto `490`;
- pérdida frente al mejor: `0.191273348886 dB`;
- ahorro de FLOP: `26.315789 %`;
- ganancia frente a GMP al mismo presupuesto: `0.382980486137 dB`.

Estos números históricos son una referencia de regresión. Para citar cifras finales debe leerse el CSV y el resumen del barrido versionado actual, no copiar números de memoria.

### 4.4 Selección automática y override explícito

`main_sweep_and_comparison.m` y `run_selected_comparison.m` usan por defecto la selección automática. Existe un override explícito para reproducir un presupuesto elegido manualmente. Cuando se usa, el informe conserva tanto el resultado automático como el valor forzado para que la decisión manual sea visible y auditable.

No se debe convertir el override en el comportamiento predeterminado ni ocultar la selección automática.

### 4.5 Figuras de artículo

La paleta IEEE utilizada es:

| Nombre | Hex | RGB |
|---|---:|---:|
| Gray | `#75787B` | `[117 120 123] / 255` |
| Blue | `#00629B` | `[0 98 155] / 255` |
| Orange | `#E87722` | `[232 119 34] / 255` |
| Green | `#00843D` | `[0 132 61] / 255` |
| Red | `#BA0C2F` | `[186 12 47] / 255` |

Convenciones visuales principales:

- GMP: línea continua y marcador circular;
- PN: línea discontinua y marcador cuadrado;
- PNNN: línea punto-raya y marcador triangular;
- referencia/original: gris y mayor grosor;
- familias ridge: color por familia y estilo distinguible por lambda;
- punto seleccionado: línea, marcador y anotación rojos;
- reducción de densidad de marcadores mediante `MarkerIndices` y equivalente TikZ `mark repeat`.

En figuras NMSE, el límite superior de Y es exactamente `-30 dB`; el inferior se redondea hacia abajo en pasos de 5 dB según los datos. La anotación identifica el punto seleccionado.

La figura de `MaxAbsRealParameter` incluye barrido principal y resultados de lambda fija. Usa escala logarítmica solo cuando todos los valores representados son positivos; si no, conserva escala lineal y lo informa.

Archivo principal de trazado:

- `toolbox/sweep/plotSweepPaperFigure.m`

### 4.6 Exportación desde una única figura MATLAB

Cada conjunto de formatos se deriva del mismo handle de figura mediante la utilidad de exportación del repositorio. Se generan:

- `.fig` para reproducción y edición MATLAB;
- `.png` a 300 dpi y fondo blanco;
- `.pdf` vectorial;
- `.tikz` mediante `matlab2tikz`.

El wrapper LaTeX independiente carga TikZ/PGFPlots/AMS math, configura compatibilidad PGFPlots 1.18 y usa variables de anchura/altura. No se reconstruye una figura diferente para cada formato.

La limpieza de archivos auxiliares LaTeX solo ocurre después de una compilación correcta. Ante fallo se conserva el log para diagnóstico.

Archivos relacionados:

- `toolbox/plotting/exportPaperFigure.m`
- `toolbox/plotting/resolveLatexmkCommand.m`
- `toolbox/plotting/runPaperFigurePreflight.m`
- `third_party/matlab2tikz/`

### 4.7 Preflight de `matlab2tikz` y LaTeX

El preflight ocurre antes de cargar datos costosos o entrenar modelos. Comprueba el submódulo, exporta una figura mínima con `matlab2tikz` y compila el wrapper.

En esta máquina, el `latexmk` de MiKTeX no podía ejecutarse directamente porque no había un Perl del sistema. La resolución implementada usa el Perl incluido con MATLAB para invocar `latexmk.pl` de MiKTeX, sin instalar ni modificar software global.

Si el preflight falla, se debe arreglar la cadena de herramientas; no se debe saltar silenciosamente para producir un conjunto parcial de formatos.

### 4.8 Espectros del punto seleccionado

Se conserva el cálculo Welch existente; no se cambió su significado matemático. Solo se generan espectros para el presupuesto seleccionado, evitando multiplicar resultados por todos los presupuestos.

La figura de error seleccionada incluye:

- espectro del objetivo;
- error de GMP;
- error de PN;
- error de PNNN.

La disposición es fija `1x2`, con leyenda compartida debajo. Los cuatro formatos de exportación salen del mismo handle.

## 5. Barrido completo versionado

El barrido validado de esquema 3 está en:

```text
results/parameter_sweep/sweep_d113e389ab78/
```

El commit `c2e4bf1` añadió 86 artefactos, entre ellos:

- `linear_sweep.mat`;
- `complexity_sweep.csv` y `.mat`;
- `fixed_lambda_linear_sweep.csv` y `.mat`;
- matrices objetivo PNNN desde `0020` hasta `0500`;
- fuente densa PNNN;
- los tres informes de selección;
- figuras principales de NMSE frente a FLOP, NMSE frente a parámetros y máximo parámetro absoluto;
- exportaciones `.fig`, `.png`, `.pdf` y `.tikz`;
- espectros del punto seleccionado en `selected_point_0340/`.

Figuras raíz versionadas:

```text
comparison_nmse_flops_sweep.{fig,png,pdf,tikz}
comparison_nmse_parameters_sweep.{fig,png,pdf,tikz}
comparison_max_abs_parameter_sweep.{fig,png,pdf,tikz}
```

En `selected_point_0340/` hay cuatro figuras espectrales, cada una en los mismos cuatro formatos.

No volver a ejecutar el barrido completo solo para comprobar que existe. Primero validar si el objetivo es reutilizarlo, regenerar únicamente presentación o hacer una corrida científica nueva.

## 6. Validaciones conocidas del estado implementado

Antes del barrido completo se registró:

- `14/14` tests rápidos aprobados;
- `checkcode` sobre 18 archivos de producción sin incidencias;
- `git diff --check` aprobado;
- preflight `matlab2tikz` + LaTeX aprobado usando el fallback de Perl de MATLAB;
- inspección visual de PNG y PDF renderizados con Poppler;
- corrección iterativa de leyendas y densidad de marcadores tras esa inspección.

Después se ejecutó y versionó el barrido completo de esquema 3 en `c2e4bf1`. Las validaciones anteriores describen la implementación estable previa y el barrido confirma que la ruta completa llegó a producir los artefactos. Cualquier cambio no confirmado posterior debe volver a probarse.

## 7. Estado no confirmado observado durante este traspaso

El working tree **no estaba limpio** porque otro chat/proceso trabajaba simultáneamente. Estos cambios no se auditaron ni se atribuyeron a esta tarea de documentación.

Se observaron modificaciones en:

```text
config/getFairDOMPComparisonConfig.m
main_sweep_and_comparison.m
run_parameter_sweep.m
run_selected_comparison.m
tests/run_main_sweep_and_comparison_quick_test.m
tests/run_operating_point_selection_tests.m
tests/run_paper_figure_integration_tests.m
tests/run_sweep_configuration_tests.m
toolbox/sweep/plotSweepPaperFigure.m
toolbox/sweep/selectOperatingPoint.m
```

También se observaron:

- eliminación aparente de `docs/paper_figure_workflow.md`;
- eliminación aparente de las rutas versionadas de los barridos históricos `sweep_0dd97cdd1cca` y `sweep_c69b734501e8`;
- nuevo árbol no versionado `results/parameter_sweep/legacy/`, posiblemente como traslado deliberado de esos barridos;
- múltiples modificaciones dentro del barrido validado `sweep_d113e389ab78`, incluidos `.mat`, figuras, CSV/TXT de selección y espectros;
- archivo nuevo no versionado `results/parameter_sweep/sweep_d113e389ab78/selected_point_0340/selected_complete_comparison.csv`.

Estas observaciones sugieren trabajo posterior de reorganización, reutilización caliente, regeneración de resultados o metadatos, pero eso es una inferencia. No asumir que las eliminaciones son accidentales ni que el traslado a `legacy/` está terminado.

Acciones prohibidas hasta coordinar ese trabajo:

- `git reset --hard`;
- `git checkout -- <archivo>` o `git restore` masivo;
- borrar `results/parameter_sweep/legacy/`;
- devolver manualmente los barridos históricos a su ubicación anterior;
- regenerar encima de `sweep_d113e389ab78`;
- hacer `git add .`, commit o push mezclando este WIP sin revisarlo;
- ejecutar formateadores o reescrituras mecánicas sobre los archivos modificados.

## 8. Cómo reanudar con seguridad

### Paso A: tomar un snapshot nuevo

```powershell
git status --short --branch
git diff --name-status
git diff --stat
git log --oneline --decorate -8
```

Comparar el resultado con las secciones 1 y 7. Si cambió, describir la diferencia antes de editar.

### Paso B: distinguir `HEAD` de WIP

Para inspeccionar la versión estable sin tocar el working tree:

```powershell
git show HEAD:config/getFairDOMPComparisonConfig.m
git show HEAD:toolbox/sweep/selectOperatingPoint.m
git show HEAD:results/parameter_sweep/sweep_d113e389ab78/operating_point_selection_summary.txt
```

Para revisar el WIP actual:

```powershell
git diff -- config/getFairDOMPComparisonConfig.m
git diff -- toolbox/sweep/selectOperatingPoint.m
git diff -- toolbox/sweep/plotSweepPaperFigure.m
git diff -- run_parameter_sweep.m run_selected_comparison.m main_sweep_and_comparison.m
git diff -- tests
```

Los `.mat`, `.fig`, PDF y PNG son binarios: revisar primero sus hashes, tamaños, fechas y los metadatos cargados por MATLAB. No concluir que cambiaron científicamente solo porque Git los marque como modificados.

### Paso C: averiguar el propósito del WIP concurrente

Confirmar con el usuario o el otro chat:

- si `legacy/` es un traslado definitivo;
- por qué se modificó el barrido ya validado;
- si `selected_complete_comparison.csv` es un entregable nuevo requerido;
- si los cambios de selección preservan el criterio científico documentado;
- si hay una ejecución MATLAB activa o recién terminada;
- qué archivos considera completos y cuáles están a medio escribir.

### Paso D: validar después de integrar

Las suites relevantes se encuentran en `tests/`. Ejecutar primero las pruebas rápidas específicas de configuración, selección, figuras e integración. Después:

1. `checkcode` sobre archivos de producción modificados;
2. `git diff --check`;
3. preflight de figuras;
4. inspección visual de PNG y PDF regenerados;
5. prueba de reutilización del barrido existente;
6. solo si es científicamente necesario, corrida completa en frío.

No aceptar únicamente “el script terminó”: comprobar resumen de selección, tablas, cuatro formatos, PDFs compilados y contenido espectral.

## 9. Mapa funcional de archivos

| Área | Archivos o directorios |
|---|---|
| Configuración | `config/getFairDOMPComparisonConfig.m` |
| Entrada completa | `main_sweep_and_comparison.m` |
| Barrido | `run_parameter_sweep.m` |
| Comparación seleccionada | `run_selected_comparison.m` |
| Firma/esquema | `toolbox/sweep/buildExperimentSignature.m` y utilidades de carga/validación |
| Selección | `toolbox/sweep/selectOperatingPoint.m` |
| Figuras de barrido | `toolbox/sweep/plotSweepPaperFigure.m` |
| Exportación | `toolbox/plotting/exportPaperFigure.m` |
| Preflight | `toolbox/plotting/runPaperFigurePreflight.m` |
| Resolución LaTeX | `toolbox/plotting/resolveLatexmkCommand.m` |
| Conversión TikZ | `third_party/matlab2tikz/` |
| Tests clave | `tests/run_sweep_configuration_tests.m`, `tests/run_operating_point_selection_tests.m`, `tests/run_paper_figure_integration_tests.m`, `tests/run_main_sweep_and_comparison_quick_test.m` |
| Barrido validado | `results/parameter_sweep/sweep_d113e389ab78/` |

## 10. Criterios de aceptación al cerrar el trabajo concurrente

Antes de declarar completada cualquier continuación, confirmar:

- rama, remoto y `HEAD` correctos;
- working tree entendido archivo por archivo;
- esquema 3 conservado o migración explícitamente justificada;
- `MaxAbsRealParameter` calculado solo sobre parámetros activos;
- advertencia de comparabilidad entre parametrizaciones conservada;
- selección automática reproducible, con tolerancia y sensibilidad reportadas;
- comparación con GMP al mismo presupuesto;
- override manual visible y no predeterminado;
- todas las figuras proceden de un único handle por conjunto de formatos;
- `.fig`, `.png`, `.pdf` y `.tikz` presentes y válidos;
- preflight anterior al trabajo costoso;
- espectros limitados al punto seleccionado;
- tests, `checkcode` y `git diff --check` aprobados;
- resultados binarios revisados, no añadidos a ciegas;
- no se perdió ningún cambio del otro chat/proceso.

## 11. Prompt listo para pegar en otro chat

```text
Trabaja exclusivamente en el repositorio:
C:\Sergi\Investigacion\Códigos\NN\PNNN_PNGMP_DOMP_fair_comparison

Lee primero, completo, CODEX_HANDOFF_2026-07-20.md. Es el documento de traspaso del trabajo de figuras de artículo, selección automática del punto operativo y MaxAbsRealParameter.

Antes de editar, informa:
- directorio actual;
- raíz Git;
- rama activa;
- commit HEAD;
- remoto;
- submódulos;
- estado completo del working tree.

El snapshot del traspaso observó la rama feature/paper-figures-selection-and-coefficient-range y HEAD c2e4bf198a94deccc9b4161c3d96dac398f7232c, pero había cambios simultáneos no confirmados. Si el estado difiere, no lo corrijas ni restaures: explica la diferencia y preserva todos los cambios.

Distingue siempre entre el estado versionado en HEAD y el WIP simultáneo descrito en la sección 7. No uses reset, checkout/restore masivo, borrado, git add ., commit ni push sin revisar el alcance. No regeneres encima del barrido validado hasta saber por qué está modificado.

Continúa la tarea que te indique usando el código, tests y resultados actuales como fuentes de verdad. Mantén intactas las decisiones científicas de split, DOMP, entrenamiento y NMSE salvo autorización explícita y justificación científica.
```

## 12. Nota sobre este propio archivo

Este Markdown fue creado únicamente como traspaso. No pretende resolver ni incorporar el WIP simultáneo. Si aparece como único archivo nuevo atribuible a este chat, puede versionarse por separado cuando el usuario decida; no debe usarse como excusa para incluir el resto del working tree en el mismo commit.
