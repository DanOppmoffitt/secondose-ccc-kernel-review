# CCC kernel alignment root-cause investigation (10x10)

Investigation-only report. No production physics path was modified.

## Question
Why is calculated dmax 48 mm while measured dmax is 12.8 mm?

## Evidence table

| Evidence | Observation | Interpretation |
|---|---|---|
| Baseline dmax | 48.00 mm | Confirms deep buildup in current model |
| Kernel-offset slope | 0.000 mm/voxel | Sensitivity to kernel depth frame |
| Convolution-offset slope | 6.000 mm/voxel | Sensitivity to indexing origin |
| Best kernel offset trial | offset=-3, dmax_diff=35.20 mm | Temporary correction potential |
| Best convolution offset trial | offset=-3, dmax_diff=17.20 mm | Temporary indexing correction potential |

## Coordinate audit summary

- Kernel depth coordinate: radial grid `r_grid_cm` starts at 0.0 and is sampled by radius.
- Kernel origin/center: interaction point at r=0; no explicit voxel-center half-step offset.
- Convolution indexing: transport loop starts at n=1, so no explicit self-dose at n=0.
- Voxel convention: world coords use `origin + index * spacing` (voxel-center convention).
- TERMA alignment: depth uses `max(d_src - SAD, 0)` in beam direction.
- Dose-grid depth indexing: CAX depth uses beam-projected source distance and keeps `>= -spacing/2`.

## Root-cause conclusion

Primary hypothesis: **kernel_model**

- If offset scans cannot close the 35.2 mm gap, the issue is likely kernel model shape (buildup model), not indexing.
- If one offset family gives near-linear expected dmax shift and approaches measured dmax, that family indicates the frame/index source.

## Recommended permanent fix

Kernel redesign or replacement with measured 6 MV kernel likely required.

Safety assessment:
- requires_kernel_redesign

## Scope guardrails

- No parameter fitting performed.
- No cohort cases run.
- No silent production-path modifications.