# simbl

`simbl` is a Python package for solving compressible boundary layer similarity equations.

## Quick Start

Install (see [Installation](installation.md) for details):

```bash
pip install similarity-bl
```

Solve a Mach 6 adiabatic flat plate:

```python
from simbl import SimilarityInputs, SolverOptions, solve_similarity

inputs = SimilarityInputs(mach_edge=6.0, temp_edge=55.0, wall_bc="adiabatic")
options = SolverOptions()

solution, info = solve_similarity(inputs, options)
print(f"f''(0) = {solution.fpp[0]:.6f}")
```

## License

BSD-3-Clause. See [LICENSE](https://github.com/uahypersonics/similarity-bl/blob/main/LICENSE) for details.
