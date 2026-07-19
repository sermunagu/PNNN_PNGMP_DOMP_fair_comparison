# Paper figures, operating-point selection, and schema v3

## External dependency

`matlab2tikz` is a Git submodule at `third_party/matlab2tikz`, pinned to
commit `806c97d99f87f8a1e99a7c54e853c25c82aac301`. After cloning, initialize it
with:

```powershell
git submodule update --init --recursive
```

MATLAB adds `third_party/matlab2tikz/src` to its path. The mandatory
preflight checks that `matlab2tikz` is callable and compiles a minimal
standalone `pgfplots` figure with `latexmk` before measurements are loaded or
training starts. On Windows/MiKTeX, if the `latexmk` launcher cannot find a
Perl script engine, the preflight and exporter use the Perl distributed with
MATLAB to invoke MiKTeX's existing `latexmk.pl`; they do not change PATH or
the system installation.

Run the preflight alone with:

```powershell
matlab -batch "addpath('config'); addpath('toolbox/sweep'); cfg=getFairDOMPComparisonConfig(pwd); disp(preflightPaperFigureToolchain(pwd,cfg.paper))"
```

## Scientific additions

Schema v3 appends `MaxAbsRealParameter` to the principal and fixed-Ridge
tables. For Complex GMP it is the maximum of the absolute real and imaginary
stored coefficient components. For PN-IQ it is the maximum across both real
coefficient vectors. For sparse PNNN it is the maximum absolute final
learnable selected by the pruning masks, including active weights and
protected biases and excluding pruned entries.

The default `near-optimal minimum-complexity criterion` uses sparse PNNN as
the reference. A point is admissible when its full-signal NMSE is within
`0.20 dB` of the best observed PNNN result and it is no worse than Complex GMP
at the same real-parameter budget. The admissible point with the fewest FLOPs
is selected, with active real parameters breaking FLOP ties. Sensitivity is
also evaluated at `[0.10 0.15 0.20 0.25] dB`.

Selection and diagnostics are written as:

- `operating_point_selection.csv`
- `operating_point_selection_summary.txt`
- `operating_point_selection_sensitivity.csv`

Selection and figure configuration are postprocessing settings and are not
part of the scientific model identity. Schema v3 itself changes the artifact
contract and therefore creates a new signed sweep directory; schema-v2
checkpoints are not migrated or silently reused.

## Figure outputs

Every public figure is exported from one figure handle to `.fig`, 300-dpi
white-background `.png`, `.tikz`, and a `.pdf` compiled from that TikZ. The
sweep root contains:

- `comparison_nmse_parameters_sweep.*`
- `comparison_nmse_flops_sweep.*`
- `comparison_max_abs_parameter_sweep.*`

Only the selected budget gets spectra, under `selected_point_NNNN`:

- `selected_output_spectrum.*`
- `selected_error_spectrum.*`
- `selected_ridge_output_spectrum.*`
- `selected_ridge_error_spectrum.*`

The standalone PDF wrapper loads `pgfplots`, `tikz`, and `amsmath`, fixes
`pgfplots` compatibility to `1.18`, defines `\figurewidth` and
`\figureheight`, and is compiled with:

```text
latexmk -pdf -interaction=nonstopmode -halt-on-error
```

Auxiliary wrapper files are removed only after success. A failed compile
preserves its `.log` and other intermediates beside the figure.

## Full workflow commands

The first schema-v3 run is cold and performs the complete `20:10:500` sweep:

```powershell
matlab -batch "run_parameter_sweep(20:10:500)"
```

After that run succeeds, the same command is the signed hot-reuse check:

```powershell
matlab -batch "run_parameter_sweep(20:10:500)"
```

The normal workflow automatically selects the operating point and generates
its four spectrum figures:

```powershell
matlab -batch "main_sweep_and_comparison"
```

An explicit manual reproduction remains available, while the workflow still
records the automatic suggestion:

```matlab
sweep = run_parameter_sweep(20:10:500);
selected = run_selected_comparison(340, sweep);
```
