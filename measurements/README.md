# Measurement contract and distribution policy

The configured article measurement is expected at the project-relative path
`measurements/experiment20260429T134032_xy.mat`.

## Required variables

The MAT-file must contain:

- `x`: finite numeric signal vector;
- `y`: finite numeric signal vector with the same number of samples as `x`.

The article data uses complex double column vectors. The loader reshapes vector
inputs to columns and interprets their roles through `mappingMode` in
`config/article_config.m`. The repository's local modeled-block X/Y convention
applies; `xy_forward` is not assigned an additional physical PA-forward meaning.

The current file also contains optional provenance fields `fs`, `description`,
and `info_signal`. They are not required by the experiment runner.

## Distribution policy

The current measurement is versioned in this research repository and is about
15.1 MB. Its external redistribution rights are not documented here. It may be
shared only with authorized tutors or project collaborators under the applicable
measurement-owner and institutional permissions. Do not publish, mirror, or
redistribute the MAT-file publicly until those rights have been confirmed.

If a distribution omits the measurement, recipients must obtain the authorized
file separately and place it at the relative path above. Do not replace it with
a different capture without documenting that change, because all reported
numerical results are measurement-specific.
