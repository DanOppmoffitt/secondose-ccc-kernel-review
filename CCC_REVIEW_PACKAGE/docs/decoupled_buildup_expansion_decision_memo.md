\# Decoupled Buildup Expansion Decision Memo



\## Status



Research-only. Candidate not frozen.



This memo summarizes the `ccc\_decoupled\_buildup\_expansion\_v1` investigation for the `TRIEXP\_DECOUPLED\_BUILDUP\_GEOMETRIC\_DILUTED\_KERNEL` research convention.



No production transport defaults were modified. No research convention was wired into the production engine router. No commissioning package was created or frozen. No patient or cohort cases were run.



\## Background



Prior CCC-native commissioning work showed that geometric dilution was necessary to recover realistic buildup behavior. Subsequent dual-exponential and tri-exponential kernel investigations improved post-dmax behavior but could not simultaneously satisfy the G1 dmax and G2 post-dmax mean gates.



The dmax sensitivity decomposition identified `longitudinal\_shape` as the dominant dmax-control lever. However, using `longitudinal\_shape` to move dmax upstream also degraded post-dmax mean agreement, indicating that the single longitudinal-shape parameter was over-coupled.



The decoupled-buildup architecture was introduced to test whether buildup and post-dmax shaping could be separated.



\## Prior Best Result



The initial decoupled-buildup probe found a partial improvement:



| Metric            |  Result |

| ----------------- | ------: |

| G1 dmax error     | 0.70 mm |

| G2 post-dmax mean |   4.19% |

| G3 post-dmax max  |   5.39% |



This preserved G1 and G3 while improving G2 relative to the prior closest longitudinal-shape compensation result, but G2 remained above the 3.0% gate.



\## Expansion Probe



The expansion probe swept around the best partial candidate.



Current best seed:



| Parameter           | Value |

| ------------------- | ----: |

| buildup\_shape       |  1.40 |

| post\_dmax\_shape     |  0.80 |

| scatter\_weight      |  0.14 |

| transition\_depth\_cm |   1.5 |

| transition\_width\_cm |   0.5 |



Expanded search axes:



| Parameter           | Values                       |

| ------------------- | ---------------------------- |

| buildup\_shape       | 1.30, 1.40, 1.50             |

| post\_dmax\_shape     | 0.80, 0.90, 1.00, 1.10, 1.20 |

| transition\_depth\_cm | 1.0, 1.5, 2.0                |

| transition\_width\_cm | 0.3, 0.5, 0.8                |

| scatter\_weight      | 0.14                         |



The run was interrupted after 110 of 135 cells. The remaining cells were in the region `buildup\_shape = 1.50` and `post\_dmax\_shape >= 1.00`, where G1 had already failed in the logged results. The outcome was therefore determined without requiring completion of the remaining cells.



\## Best Observed Candidates



| Rank | buildup\_shape | post\_dmax\_shape | transition\_depth\_cm | transition\_width\_cm | G1 dmax error | G2 mean | G3 max | G1   | G2   | G3   |

| ---: | ------------: | --------------: | ------------------: | ------------------: | ------------: | ------: | -----: | ---- | ---- | ---- |

|    1 |          1.50 |            0.80 |                 1.5 |                 0.3 |       0.70 mm |   4.06% |  5.26% | PASS | FAIL | PASS |

|    2 |          1.40 |            0.80 |                 1.5 |                 0.3 |       0.70 mm |   4.10% |  5.31% | PASS | FAIL | PASS |

|    3 |          1.50 |            0.80 |                 1.5 |                 0.5 |       0.70 mm |   4.16% |  5.36% | PASS | FAIL | PASS |



Best observed G2 post-dmax mean was 4.06%, still above the 3.0% gate.



\## Findings



The decoupled-buildup architecture improves the G2 error relative to the prior closest longitudinal-shape compensation result, but the improvement plateaus near 4.0%.



Increasing `post\_dmax\_shape` above 0.80 did not close G2. Instead, it consistently moved the solution toward G1 failure by shifting dmax away from the measured 12.8 mm target.



The best result was boundary-pinned at:



\* upper tested `buildup\_shape`

\* lower tested `post\_dmax\_shape`

\* lower tested `transition\_width\_cm`



However, the direction implied by this boundary pinning does not support a simple extension of the same architecture. The low `post\_dmax\_shape` side preserves G2 better but fails to move dmax sufficiently unless paired with stronger buildup shaping. Increasing post-dmax shaping moves dmax unfavorably and breaks G1.



Transition-depth and transition-width tuning produced only modest G2 changes and did not approach the 3.0% G2 gate.



\## Decision



`TRIEXP\_DECOUPLED\_BUILDUP\_GEOMETRIC\_DILUTED\_KERNEL` should not be frozen.



Decision: \*\*FAILURE / DO NOT FREEZE\*\*



The architecture produced a useful partial improvement but did not achieve simultaneous G1, G2, and G3 agreement within the tested research space.



\## Interpretation



The investigation shows that buildup and post-dmax behavior are partially separable, but not sufficiently separable using the current decoupled-shape formulation.



The remaining G2 gap appears structurally limited near approximately 4.0% when G1 and G3 are preserved.



This result argues against further broad expansion of the same decoupled-buildup parameter space.



\## Recommended Next Direction



Do not continue expanding ordinary buildup-shape, post-dmax-shape, transition-depth, or transition-width sweeps.



The next research step should investigate a different correction strategy that preserves the recovered dmax behavior while correcting the residual post-dmax mean error.



A reasonable next hypothesis is a tightly bounded research-only post-dmax residual correction applied after buildup behavior is established, with explicit preservation of:



\* G1 dmax placement

\* G3 max-error behavior

\* 10 cm absolute calibration anchor

\* production-path isolation



The goal should be to reduce post-dmax mean error from approximately 4.06% to ≤3.0% without moving dmax.



\## Production Isolation



This work remains research-only and candidate\_not\_frozen.



No production defaults were modified. The production engine router remains limited to:



```python

\["analytical", "ccc"]

```



No clinical validation, commissioning, or production-readiness claim is made from this probe.



