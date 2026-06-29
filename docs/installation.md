# Installation

## From PyPI

```bash
pip install similarity-bl
```

To upgrade an existing installation:

```bash
pip install --upgrade similarity-bl
```

## From Source

```bash
git clone https://github.com/uahypersonics/similarity-bl.git
cd similarity-bl
pip install -e .
```

## Optional Extras

For the `simbl` CLI:

```bash
pip install "similarity-bl[cli]"
```

Includes:

- [typer](https://typer.tiangolo.com/) for the command-line interface

For development (tests, linting, and docs):

```bash
pip install -e ".[dev]"
```

Includes:

- [pytest](https://docs.pytest.org/) and [pytest-cov](https://pytest-cov.readthedocs.io/) for testing
- [ruff](https://docs.astral.sh/ruff/) for linting
- [zensical](https://zensical.org/) for building the documentation

## Verify Installation

```bash
simbl --version
```

Or in Python:

```python
import simbl
print(simbl.__version__)
```
