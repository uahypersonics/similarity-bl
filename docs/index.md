# simbl

`simbl` is a Python package for solving compressible boundary layer similarity equations.

## Quick Start

### Install

```bash
pip install similarity-bl
```

### Run

=== "CLI"

    ```bash
    simbl solve --mach 6.0 --temp-edge 55.0 --wall adiabatic
    ```

=== "API"

    ```python
    from simbl import SimilarityInputs, SolverOptions, solve_similarity

    inputs = SimilarityInputs(mach_edge=6.0, temp_edge=55.0, wall_bc="adiabatic")
    options = SolverOptions()

    solution, info = solve_similarity(inputs, options)
    print(f"f''(0) = {solution.fpp[0]:.6f}")
    ```

## Feedback & Contributing

Questions, bug reports, and contributions are welcome. If something unexpected
comes up while using this solver, or there are ideas for improvement, opening
an issue or starting a discussion is the best first step.

Using a label when opening an issue helps prioritize and track requests:

- [Ask a question](https://github.com/uahypersonics/similarity-bl/issues/new?labels=question)
- [Report a bug](https://github.com/uahypersonics/similarity-bl/issues/new?labels=bug)
- [Suggest a feature](https://github.com/uahypersonics/similarity-bl/issues/new?labels=enhancement)

For those interested in contributing code or documentation, the
[Contributing Guide](https://github.com/uahypersonics/similarity-bl/blob/main/CONTRIBUTING.md)
covers how to set up a dev environment, run tests, and submit a PR.

## License

BSD-3-Clause. See [LICENSE](https://github.com/uahypersonics/similarity-bl/blob/main/LICENSE) for details.
