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

## Map Onto A Structured Grid

Write the similarity solution as JSON, then generate the baseflow mapping
configuration:

```bash
simbl solve simbl_config.toml --output simbl.json
simbl baseflow init
```

Edit `simbl_baseflow.toml` to identify the solution, edge state, input grid, and
output file:

```toml
[baseflow]
solution_input = "simbl.json"
edge_state_input = "edge_state.json"
grid_input = "grid.hdf5"
output = "baseflow.hdf5"
similarity_transform = "levy_lees"
geometry_transform = "none"
farfield_fill = "edge"
endpoint_tolerance = 1.0e-3
dtype = "f8"
```

All artifact paths are resolved relative to the configuration file. Run the
mapping with:

```bash
simbl baseflow run simbl_baseflow.toml
```

The command preserves the input coordinates and writes `uvel`, `vvel`, `wvel`,
`temp`, `pres`, `dens`, and `visc` fields through `cfd-io`. Cubic splines map the
solution inside its eta domain. Grid points beyond the mapped `eta_max` receive
the exact freestream values from the edge FlowState rather than extrapolated
spline values.

Farfield filling is allowed only when the similarity endpoint satisfies the
configured convergence tolerance for `fp`, `tau`, `fpp`, and `taup`. The current
workflow supports planar Falkner-Skan profiles. Mangler geometry mapping and
Falkner-Skan-Cooke crossflow profiles are rejected explicitly until those
transformations are implemented.

See the [Similarity to Physical Coordinate Transform](../theory/similarity_to_physical_coordinate_transform/index.md)
theory page for definitions.
