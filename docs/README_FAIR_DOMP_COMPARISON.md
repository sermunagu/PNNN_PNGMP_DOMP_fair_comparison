# Fair PNNN versus PN-GMP DOMP comparison

This project compares Complex GMP, coupled and independent PN-IQ models, and
PNNN under one modeled-block X/Y convention. Support selection, regularization,
and neural schedules use only the identification domain. Final metrics are
reported over identification and over the complete signal.

## Main protocol

- Identification uses a deterministic 10% of the capture (seed 1004).
- Internal train/validation partitions identification 85/15 (seed 42).
- DOMP paths and neural checkpoints use internal train/validation only.
- Final coefficients and neural refits use all identification samples.
- Full-signal rows are used only for final evaluation.

Full-signal evaluation includes the identification samples. It is not an
independent hold-out test and must not be interpreted as a generalization set.
The experiment signature hashes the pre-DC X/Y arrays and every reuse-critical
configuration field before any saved selection can be reused.

## Linear and neural models

The normal comparison retains three distinct linear estimators on fixed
supports:

1. LS (`lambda=0`);
2. fixed Ridge (`lambda=1e-5`);
3. Ridge selected from the configured grid using internal validation only.

Fixed Ridge changes coefficient estimation, not support, parameter count, or
inference FLOPs. In the completed full-signal run for this capture, fixed Ridge
did not improve the two 344-parameter main linear models: Independent PN-IQ
changed from -38.839 dB (LS) to -37.120 dB, and parameter-matched Complex GMP
from -38.458 dB to -36.860 dB. Internal validation selected `lambda=0` for both.

The PNNN uses periodic memory `M=13`, orders `[1 3 5 7]`, full features,
sigmoid activation, and two real outputs. Its input dimension is 84. The N12
dense source has 1046 trainable real scalars. Sparse targets are produced by
independent global magnitude pruning from that same dense source: biases are
protected, zero weights remain frozen, and every target reuses the same seed
and historical fine-tuning budget.

For the historical 344-parameter sparse N12 point, the ideal sparse count is
848 FLOPs/sample when zero weights are skipped; executing the original dense
matrices costs 2252 FLOPs/sample. Both are analytical operation counts, not
MATLAB timings or hardware measurements. Sigmoid evaluations, exponentials,
square roots, divisions, and magnitudes remain separately reported operations.

## Parameter sweep

`run_parameter_sweep` compares three genuine curves on one parameter grid:

- Complex GMP uses one maximum DOMP path on internal train and one on all
  identification rows; every point is an exact prefix.
- Independent PN-IQ uses its own genuine PN-DOMP paths and exact `P/2` feature
  prefixes after one structural reduction per domain.
- Sparse PNNN N12 prunes every target independently from one signed dense N12
  source; pruning is never progressive.

The historical 344-parameter Independent PN-IQ point is also reported, but it
is a separate marker. Its features are induced by the Complex GMP DOMP-100
support and are not forced into the genuine PN-DOMP path. A difference between
the two 344-parameter PN-IQ points therefore reflects different feature
selection protocols, not a regression.

Sweep outputs are written below:

```text
results/parameter_sweep/sweep_<identity-digest>/
    linear_sweep.mat
    sweep_dense_source.mat
    pnnn_target_0150.mat
    ...
    complexity_sweep.csv
    complexity_sweep.mat
    comparison_nmse_parameters_sweep.png
```

The linear family is one atomic resume unit. Each PNNN target is one atomic,
deterministically named artifact, so a corrupt or incompatible target causes
only that target to be recalculated. The consolidated MAT contains tables,
signatures, artifact names, summarized linear supports, and the 344 comparison;
it does not duplicate all PNNN predictions.

## Commands

From the project root:

```powershell
matlab -batch "run('tests/run_domp_unit_test.m')"
matlab -batch "run('tests/run_experiment_signature_test.m')"
matlab -batch "run('tests/run_linear_complexity_sweep_test.m')"
matlab -batch "run('tests/run_pnnn_shared_dense_sweep_test.m')"
matlab -batch "run('tests/run_fair_comparison_smoke_test.m')"
matlab -batch "run_fair_PNNN_vs_PNGMP_DOMP"
matlab -batch "run_parameter_sweep"
```

The sweep is disabled in the normal comparison configuration. The full sweep
and the normal full-measurement pipeline are intentionally separate entry
points. Neither analytical FLOPs nor sparse-kernel assumptions establish
latency, energy, FPGA utilization, or memory-bandwidth superiority.
