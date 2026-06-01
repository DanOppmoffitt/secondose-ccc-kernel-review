# Longitudinal-shape compensation probe (research-only)

**Status:** candidate_not_frozen / research_only. Production NOT modified. primary_decay bound NOT relaxed.

- Probe: `ccc_longitudinal_compensation`
- Measured dmax: 12.8 mm
- Gates: G1 <= 2.0 mm, G2 <= 3.0 %, G3 <= 8.0 %

## Decision

G1_RECOVERED_BUT_G2G3_SACRIFICED — longitudinal_shape in [1.1, 1.5] recovers G1 but no scatter_weight level in the probed band pulls G2/G3 fully under gate. Closest cell: longitudinal_shape=1.10 + scatter_weight=0.30 (G2=5.608%, G3=6.800%). scatter_weight compensation is insufficient on its own; an additional orthogonal lever is required. Candidate NOT frozen.

## Compensation grid (longitudinal_shape x scatter_weight)

| long | scatter | dmax mm | G1 err mm | G1 | G2 mean % | G2 | G3 max % | G3 | all |
|------|---------|---------|-----------|----|-----------|----|----------|----|-----|
| 1.10 | 0.14 | 13.5 | 0.7 | ✓ | 5.7467 | ✗ | 6.925 | ✓ | ✗ |
| 1.10 | 0.22 | 13.5 | 0.7 | ✓ | 5.6761 | ✗ | 6.8615 | ✓ | ✗ |
| 1.10 | 0.30 | 13.5 | 0.7 | ✓ | 5.6077 | ✗ | 6.8 | ✓ | ✗ |
| 1.10 | 0.38 | 15.0 | 2.2 | ✗ | 5.5428 | ✗ | 6.7418 | ✓ | ✗ |
| 1.20 | 0.14 | 13.5 | 0.7 | ✓ | 6.062 | ✗ | 7.2251 | ✓ | ✗ |
| 1.20 | 0.22 | 13.5 | 0.7 | ✓ | 5.9988 | ✗ | 7.1682 | ✓ | ✗ |
| 1.20 | 0.30 | 13.5 | 0.7 | ✓ | 5.9376 | ✗ | 7.1131 | ✓ | ✗ |
| 1.20 | 0.38 | 13.5 | 0.7 | ✓ | 5.8783 | ✗ | 7.0596 | ✓ | ✗ |
| 1.30 | 0.14 | 13.5 | 0.7 | ✓ | 6.3306 | ✗ | 7.5076 | ✓ | ✗ |
| 1.30 | 0.22 | 13.5 | 0.7 | ✓ | 6.274 | ✗ | 7.4499 | ✓ | ✗ |
| 1.30 | 0.30 | 13.5 | 0.7 | ✓ | 6.2193 | ✗ | 7.394 | ✓ | ✗ |
| 1.30 | 0.38 | 13.5 | 0.7 | ✓ | 6.1662 | ✗ | 7.3398 | ✓ | ✗ |
| 1.40 | 0.14 | 13.5 | 0.7 | ✓ | 6.5616 | ✗ | 7.7491 | ✓ | ✗ |
| 1.40 | 0.22 | 13.5 | 0.7 | ✓ | 6.511 | ✗ | 7.6975 | ✓ | ✗ |
| 1.40 | 0.30 | 13.5 | 0.7 | ✓ | 6.462 | ✗ | 7.6475 | ✓ | ✗ |
| 1.40 | 0.38 | 13.5 | 0.7 | ✓ | 6.4145 | ✗ | 7.599 | ✓ | ✗ |
| 1.50 | 0.14 | 13.5 | 0.7 | ✓ | 6.7619 | ✗ | 7.9573 | ✓ | ✗ |
| 1.50 | 0.22 | 13.5 | 0.7 | ✓ | 6.7166 | ✗ | 7.9111 | ✓ | ✗ |
| 1.50 | 0.30 | 13.5 | 0.7 | ✓ | 6.6727 | ✗ | 7.8663 | ✓ | ✗ |
| 1.50 | 0.38 | 13.5 | 0.7 | ✓ | 6.6302 | ✗ | 7.8229 | ✓ | ✗ |

_Research-only. ccc_longitudinal_compensation probe, candidate_not_frozen. No production integration, no router changes, no freeze, no patient/cohort run, no validation claim. Production-adjacent primary_decay bound NOT relaxed. TRIEXP_GEOMETRIC_DILUTED_KERNEL base held fixed; two-axis (longitudinal_shape x scatter_weight) compensation grid only._
