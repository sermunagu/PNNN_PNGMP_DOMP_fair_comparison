# Full-signal PNNN versus PN-GMP DOMP comparison

This study compares complex GMP, phase-normalized GMP formulations, and PNNN
under the identification protocol used by the original GMP work. The validated
main run is:

```text
results/full_signal_domp_comparison/20260715_142905/
```

The earlier disjoint-test experiment remains unchanged under
`results/fair_domp_comparison/20260714_013938/` as an exploratory historical
control; it is not the main result.

## Protocol

- `identificationIndices`: the deterministic original 10% selector, seed 1004,
  49,151 samples;
- `fullSignalIndices`: all 491,520 samples, including every identification
  sample;
- `internalTrainIndices` and `internalValidationIndices`: an internal split of
  identification used only to select lambda, epochs, and fine-tuning duration.

After hyperparameter selection, every final linear model regenerates its DOMP
support and fits its coefficients on all 49,151 identification samples. Every
final PNNN is restarted from the same seed, computes normalization statistics on
all identification samples, and is trained for the fixed selected epoch count on
those same samples. The complete signal is used only for final evaluation.

**The full-signal evaluation includes the identification samples and therefore
must not be interpreted as an independent generalization test.**

All models use the same complex temporal NMSE:

```text
10*log10(sum(abs(y-yhat).^2)/sum(abs(y).^2))
```

## Models and parameter matching

DOMP is the only support selector. It runs during identification and is not part
of inference. The coupled PN-IQ formulation reuses the complex GMP DOMP-100
support and is numerically equivalent to it; the measured relative prediction
error is `1.325559e-11`.

The principal comparison fixes exactly 344 active real parameters:

| Model | Real parameters | Identification NMSE (dB) | Full-signal NMSE (dB) | FLOPs/sample |
|---|---:|---:|---:|---:|
| Independent PN-IQ | 344 | -39.004685 | -38.838825 | 1026 |
| Complex GMP parameter-matched | 344 | -38.672161 | -38.458395 | 1837 |
| PNNN N12 sparse | 344 | -35.696928 | -35.659914 | 848 |

Independent PN-IQ has the lowest NMSE in this capture. It improves the
parameter-matched GMP by 0.380 dB and the sparse PNNN by 3.179 dB on the full
signal. The sparse PNNN has the lowest arithmetic count, provided an
implementation actually skips its exactly zero weights. This is a trade-off,
not a claim of universal or hardware superiority.

Secondary controls remain in `comparison_results.csv`: complex GMP DOMP-100,
coupled PN-IQ, independent I/Q without phase normalization, reduced independent
PN-IQ, PNNN H4 dense, and PNNN N12 dense.

## Final PNNN refit

The PNNN input has `D=84` phase-normalized features, one sigmoid hidden layer, and
two real outputs. A width-`H` network contains

```text
H*(D + 3) + 2
```

real scalars, including weights and biases. Internal selection chose epoch 1259
for H4 and epoch 1258 for N12. The final N12 sparse model starts from the final
identification-trained N12, protects all 14 biases, retains 330 weights, freezes
the zero mask, and fine-tunes for 163 fixed epochs on all identification
samples. It therefore has 344 active parameters and 68.0233% weight sparsity.

## FLOPs/sample

`FLOPs/sample` counts real additions and multiplications required to produce one
complex output sample, including used-feature generation, coefficient/weight
application, accumulation, biases, phase normalization, and phase restoration.
The convention is one FLOP per real addition or multiplication, two per complex
addition, and six per complex multiplication.

Operations without an agreed conversion are reported separately. For example,
Independent PN-IQ also requires 12 magnitudes (including 12 square roots) and 2
divisions; N12 sparse requires 14 magnitudes, 2 divisions, and 12 sigmoid
activations (up to 12 exponentials), in addition to phase normalization and
restoration arithmetic already included in its FLOP count.

The sparse-PNNN figure of 876 FLOPs/sample assumes zero weights are skipped. An
implementation that executes the original full matrices requires 2252
FLOPs/sample. FLOPs do not directly determine MATLAB runtime, FPGA resources,
latency, power, or energy.

## Commands

From the project root:

```powershell
matlab -batch "run('tests/run_fair_comparison_smoke_test.m')"
matlab -batch "run('run_fair_PNNN_vs_PNGMP_DOMP.m')"
```

The smoke test must finish with:

```text
FULL-SIGNAL DOMP COMPARISON SMOKE TEST: PASS
```

The complete run creates a timestamped directory below
`results/full_signal_domp_comparison/` containing the CSV/MAT summaries, DOMP
supports, full-signal predictions, figures, metadata, and the two-page English
report. It never overwrites the historical disjoint-test result.

## Validated artifacts

The main numerical sources are:

- `comparison_results.csv` and `comparison_results.mat`;
- `parameter_summary.csv`;
- `complexity_flops.csv`;
- `selected_hyperparameters.csv`;
- `split_indices.mat` and `comparison_config.mat`;
- `domp_supports.mat`;
- `full_signal_predictions.mat`.

The long Spanish guide is `docs/PN_GMP_PNNN_defense_guide.pdf`; the short English
report is stored inside the validated result directory as
`fair_comparison_report.pdf`.
