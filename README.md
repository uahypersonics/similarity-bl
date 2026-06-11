# simbl

Compressible similarity boundary layer solver with Falkner-Skan and Falkner-Skan-Cooke support.

[![Test](https://github.com/uahypersonics/similarity-bl/actions/workflows/test.yml/badge.svg)](https://github.com/uahypersonics/similarity-bl/actions/workflows/test.yml)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20648457.svg)](https://doi.org/10.5281/zenodo.20648457)
[![PyPI](https://img.shields.io/pypi/v/similarity-bl)](https://pypi.org/project/similarity-bl/)
[![Docs](https://img.shields.io/badge/docs-mkdocs-blue)](https://uahypersonics.github.io/similarity-bl/)
[![License](https://img.shields.io/badge/License-BSD_3--Clause-blue.svg)](LICENSE)
[![Python](https://img.shields.io/badge/python-≥3.11-blue.svg)](https://www.python.org/downloads/)
[![Ruff](https://img.shields.io/endpoint?url=https://raw.githubusercontent.com/astral-sh/ruff/main/assets/badge/v2.json)](https://github.com/astral-sh/ruff)

## Install

```bash
pip install similarity-bl
```

## Quick Start

**API:**

```python
from simbl import SimilarityInputs, SolverOptions, solve_similarity

inputs = SimilarityInputs(mach_edge=6.0, temp_edge=55.0, wall_bc="adiabatic")
options = SolverOptions()

solution, info = solve_similarity(inputs, options)
print(f"f''(0) = {solution.fpp[0]:.6f}")
print(f"g(0)   = {solution.g[0]:.6f}")
```

**CLI:**

```bash
simbl init                   # generate a template config file
simbl solve config.toml      # solve from config
simbl solve --mach 6.0 --temp-edge 55.0 --wall adiabatic
```

## Documentation

Full documentation: https://uahypersonics.github.io/similarity-bl

## Citation

If you use `simbl` in your research, please cite it:

```bibtex
@software{simbl,
  author = {Hader, Christoph},
  title = {simbl: Compressible Similarity Boundary Layer Solver},
  url = {https://github.com/uahypersonics/similarity-bl},
  year = {2026}
}
```

## Code Style

This project follows established Python community conventions so that
contributors can focus on the physics rather than inventing formatting rules.

| Convention | What it covers | Reference |
|---|---|---|
| [PEP 8](https://peps.python.org/pep-0008/) | Code formatting, naming, whitespace | Python standard style guide |
| [PEP 257](https://peps.python.org/pep-0257/) | Docstring structure (triple-quoted, imperative mood) | Python standard docstring conventions |
| [numpydoc](https://numpydoc.readthedocs.io/en/latest/format.html) | Docstring sections (`Parameters`, `Returns`, `Attributes`) | NumPy/SciPy docstring standard — the norm for scientific Python |
| [Ruff](https://docs.astral.sh/ruff/) | Automated linting and formatting | Enforces PEP 8 compliance automatically |
| [typing / TYPE_CHECKING](https://docs.python.org/3/library/typing.html#typing.TYPE_CHECKING) | Type hints for IDE support and static analysis | Python standard library |

## Releasing

This project uses [Semantic Versioning](https://semver.org/) (`vMAJOR.MINOR.PATCH`):

- **MAJOR** (`v1.0.0`, `v2.0.0`): Breaking API changes
- **MINOR** (`v0.3.0`, `v0.4.0`): New features, backward-compatible
- **PATCH** (`v0.3.1`, `v0.3.2`): Bug fixes, minor corrections

To publish a new version to [PyPI](https://pypi.org/project/similarity-bl/):

1. Commit and push to `main`
2. Tag and push:
   ```bash
   git tag -a vMAJOR.MINOR.PATCH -m "Release vMAJOR.MINOR.PATCH"
   git push origin vMAJOR.MINOR.PATCH
   ```
