\# SeconDose CCC External Review Request



\## Objective



This package is for external technical review of the research-only CCC-native commissioning work in SeconDose.



Please review the physics flow and commissioning evidence. The main question is:



\*\*Why do multiple CCC kernel/shape architectures converge near a persistent \~4% post-dmax mean residual while preserving acceptable dmax and G3 behavior?\*\*



\## Current best candidate behavior



Measured TrueBeam 6 MV 10x10 water PDD:



\- measured dmax = 12.8 mm



Best current research candidate:



\- G1 dmax error ≈ 0.7 mm, PASS

\- G2 post-dmax mean ≈ 4.0%, FAIL

\- G3 post-dmax max ≈ 5.3%, PASS



\## Gates



\- G1: dmax error ≤ 2.0 mm

\- G2: post-dmax mean error ≤ 3.0%

\- G3: post-dmax max error ≤ 8.0%



\## What has already been investigated



\- Missing geometric dilution: confirmed and corrected in research path

\- Dual-exponential kernel family

\- Tri-exponential kernel family

\- Sub-2.0 primary decay probe

\- 1.5 mm resolution confirmation

\- Proximal shift correction

\- dmax sensitivity decomposition

\- longitudinal\_shape compensation

\- decoupled buildup/post-dmax shaping

\- post-dmax residual correction



\## Current interpretation



Geometric dilution was necessary and produced the largest improvement.



Longitudinal shape controls dmax, but coupling to post-dmax behavior prevents simultaneous G1 and G2 closure.



Decoupled buildup/post-dmax shaping improved G2 but plateaued around 4%.



Post-dmax residual correction preserved dmax but produced only minimal additional G2 improvement.



\## Review questions



1\. Does the persistent \~4% post-dmax mean residual suggest a deficiency in:

&#x20;  - TERMA generation?

&#x20;  - scatter generation?

&#x20;  - cone transport?

&#x20;  - aperture / fluence modeling?

&#x20;  - normalization / calibration?

&#x20;  - measured-data interpretation?

&#x20;  - something else?



2\. Would you continue kernel-family development, or redirect to a different CCC subsystem?



3\. Is there a known CCC failure mode that produces:

&#x20;  - realistic dmax,

&#x20;  - acceptable maximum post-dmax error,

&#x20;  - but persistent mean post-dmax bias?



4\. What specific code path or physics assumption should be reviewed next?

