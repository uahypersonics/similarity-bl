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
