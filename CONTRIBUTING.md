# Contributing to simbl

Contributions are welcome, whether in the form of bug fixes, documentation improvements,
or new physics. This guide covers everything needed to get started.

## Setting Up a Development Environment

Clone the repository and install in editable mode with all development dependencies:

```
git clone https://github.com/uahypersonics/similarity-bl.git
cd similarity-bl
pip install -e ".[dev,cli]"
```

This installs the package in editable mode along with the test suite, linter, and
documentation tools listed under `[project.optional-dependencies]` in `pyproject.toml`.

## Running the Tests

```
pytest
```

For a coverage report:

```
pytest --cov=simbl --cov-report=term-missing
```

All tests must pass before submitting a pull request.

## Code Style

This project uses [Ruff](https://docs.astral.sh/ruff/) for linting and formatting.
Run it before committing:

```
ruff check .
ruff format .
```

The CI pipeline runs both checks automatically. See the code style table in
[README.md](README.md) for the full conventions (PEP 8, PEP 257, numpydoc).

## Building the Documentation

```
mkdocs serve
```

This starts a local server at `http://127.0.0.1:8000` where you can preview changes
before pushing. The documentation is built and deployed automatically on merge to `main`.

## Submitting Changes

1. Fork the repository and create a branch from `main`:

```
git checkout -b branch-name
```

2. Make changes, add tests if applicable, and ensure `pytest` and `ruff` both pass.

3. Open a pull request against `main` with a clear description of what the change does
and why.

One feature or fix per pull request. Keep PRs focused so they are easy to review.

## Contributing a VnV Case

Validation cases are provided in `vnv/`. To add a new one:

1. Create a directory using the naming convention:
   `<geometry>_mach_<M>_re1_<Re>_<wall_bc>[_tw_<Tw>]`
2. Add the following structure:
   ```
   data/        reference baseflow or tabulated data
   figures/     pre-generated comparison figures
   scripts/     validate.py - reproduces all figures
   README.md    case description, conditions, expected outcome
   ```
3. Add a row to the cases table in [`vnv/README.md`](vnv/README.md).

See an existing case (e.g. `vnv/sharp_flat_plate_mach_05pt00_re1_11pt40e6_isothermal_tw_300_k/`)
for a concrete reference.

## Contributing Physics or Equations

Changes that touch the governing equations, transformations, or numerical methods should:

- Reference the relevant literature or other sources as
  appropriate

- Be consistent with the derivations documented at
  [uahypersonics.github.io/similarity-bl](https://uahypersonics.github.io/similarity-bl)

- Include a validation case or test that demonstrates correctness against a known result

If a proposed change may be inconsistent with the formulation, open an issue first to
discuss before writing code.

## Reporting Bugs

Open an issue on GitHub with:

- A minimal reproducible example

- The `simbl` version (`pip show similarity-bl`)

- Python version and operating system

- Expected vs. actual behavior

## License

All contributions are made under the [BSD-3-Clause License](LICENSE).
