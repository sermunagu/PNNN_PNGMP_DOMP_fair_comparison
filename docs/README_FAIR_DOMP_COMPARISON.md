# Fair PNNN versus PN-GMP DOMP comparison

This project compares PNNN with the already validated PN-GMP/DOMP study on
one frozen partition. The linear results, polynomial supports, metrics, and
TEST predictions are loaded from the validated snapshot; the current runner
does not execute DOMP or refit any linear model.

The common domains are:

- `FIT_POOL`: deterministic classic 4% selection, seed 1004;
- `TEST`: exact disjoint complement of `FIT_POOL`;
- `TRAIN`: amplitude-stratified 85% of `FIT_POOL`, seed 42;
- `VAL`: remaining 15% of `FIT_POOL`.

For the new PNNN path, normalization and gradient updates use only `TRAIN`.
The best dense and fine-tuned checkpoints are selected only through `VAL`.
`FIT_POOL` is reported as the union of `TRAIN` and `VAL`; it is not used for a
second PNNN refit. `TEST` is evaluated once after each final PNNN model has
been frozen. The saved split is loaded directly from the validated linear
snapshot so PNNN cannot create a private partition.

## PNNN models

The phase-normalized PNNN input is built with periodic memory `M=13`, orders
`[1 3 5 7]`, full features, sigmoid activation, and two real outputs. The active builder
currently gives `D=84`. A one-hidden-layer network with width `H` has

```text
H*(D + 3) + 2
```

trainable real scalars, including every weight and bias. Three models are
trained from scratch with seed 42 and no warm start:

- `PNNN H4 dense`: 350 parameters, the small dense control;
- `PNNN N12 dense`: 1046 parameters before pruning;
- `PNNN N12 sparse`: N12 globally pruned by weight magnitude to the dynamic
  `NumRealParameters` of `Independent PN-IQ full` (358 in the validated run).

N12 pruning is global and unstructured. Biases are protected, masked weights
are frozen, and fine-tuning uses `TRAIN` with `VAL` checkpoint selection. The
runner asserts the exact active count and exact zero-mask integrity after
fine-tuning. The historical N25 result is not retrained by default.

The main parameter comparison is:

- `Independent PN-IQ full`;
- `Complex GMP DOMP parameter-matched`;
- `PNNN N12 sparse`.

A wider sparse network may retain a more useful representation than a narrow
dense network with a similar active-parameter count. This is an experimental
hypothesis tested by H4 versus sparse N12, not a general guarantee.

## Update-scaled training budget

Keeping 150 epochs would sharply reduce optimization updates because the fair
`TRAIN` set is much smaller than the historical 70% PNNN training set. The
runner therefore converts the historical schedule to update counts and scales
maximum epochs, learning-rate drop period, validation patience, and pruning
fine-tuning epochs. Batch size remains 1024. Results record:

- `IterationsPerEpoch`;
- `DenseTrainingUpdates`;
- `FineTuneUpdates`;
- `BestDenseEpoch`;
- `BestFineTuneEpoch`.

No TEST sample participates in normalization, training, pruning, fine-tuning,
early stopping, or checkpoint selection.

## NMSE and sparse FLOPs

Every model retains the common complex temporal NMSE definition:

```text
10*log10(sum(abs(y-yhat).^2)/sum(abs(y).^2)).
```

For `PNNN N12 sparse`, two core arithmetic counts are reported:

- `DenseExecutionCoreFLOPsPerSample`: full N12 matrix execution, including
  multiplications by stored zeros;
- `IdealSparseCoreFLOPsPerSample`: only nonzero weights, requiring a sparse
  kernel that actually skips zero products.

The ideal count is not a guaranteed MATLAB or hardware cost. Activation evaluations, worst-case
exponential, square root, division, and absolute-value calls remain separate
special-operation columns. These algorithmic counts do not directly establish
FPGA resources, latency, power, or wall-clock runtime.

## Commands

From the project root, run the smoke test first and then the fair PNNN update:

```powershell
matlab -batch "run('tests/run_fair_comparison_smoke_test.m')"
matlab -batch "run('run_fair_PNNN_vs_PNGMP_DOMP.m')"
```

The full run creates one timestamped directory below
`results/fair_domp_comparison/`. It writes updated comparison CSVs, both plots,
the LaTeX/PDF report, metadata, and TEST-only predictions. The frozen linear
snapshot remains unchanged.

## Latest completed PNNN run

The PNNN-only update completed in
`results/fair_domp_comparison/20260714_013938/`. The validated linear rows were
loaded from `20260714_000828` and were not recomputed.

| Model | Active params | TEST NMSE (dB) | Dense core FLOPs | Ideal sparse core FLOPs |
|---|---:|---:|---:|---:|
| PNNN H4 dense | 350 | -33.459179 | 876 | 876 |
| PNNN N12 dense | 1046 | -36.469873 | 2252 | 2252 |
| PNNN N12 sparse | 358 | -35.457309 | 2252 | 876 |

N12 sparse reached the dynamic target exactly: 344 active weights plus 14
protected biases, for 358 active parameters and 66.6667% weight sparsity. H4
and N12 dense each used 50,400 gradient updates. Sparse fine-tuning used 6,732
updates; its best VAL checkpoint was epoch 387. The ideal sparse value remains
an analytical count and requires an execution kernel that skips zero weights.
