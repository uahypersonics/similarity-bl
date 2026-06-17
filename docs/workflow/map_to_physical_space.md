# Optionally Map To Physical Space

Similarity profiles are naturally functions of $\eta$. For comparison against
CFD, experiment, or geometry-based data, map the profile to the physical
wall-normal coordinate $y$.

The inverse transform uses the temperature profile $\tau = T/T_e$ and local edge
scales.

| Equation family | Default inverse transform |
|---|---|
| `falkner_skan` | Levy-Lees |
| `falkner_skan_cooke` | Illingworth-Stewartson |

Physical-space mapping is optional. If the downstream task only needs similarity
profiles, the eta-space solution can be used directly.

See the [Similarity to Physical Coordinate Transform](../theory/similarity_to_physical_coordinate_transform/index.md)
theory page for definitions.
