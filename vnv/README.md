# Validation and Verification (VnV)

Reproducible cases that confirm `simbl` is correct.

- **Validation** — compares the similarity solution against independent reference data
  (CFD, experiment, published tables). Answers: does the model represent reality?
- **Verification** — checks the numerical implementation against exact or analytic results.
  Answers: is the code solving the equations correctly?

Each case is self-contained: it includes the reference data, a script that
reproduces all figures, and pre-generated figures so the outcome is visible
without running anything.

## Verification

*(cases coming soon)*

## Validation

| Case | Geometry | Mach | Re1 (1/m) | Wall BC | Solver | Reference |
|---|---|---|---|---|---|
| [sharp_flat_plate_mach_05pt00_re1_11pt40e6_isothermal_tw_300_k](validation/falkner_skan/sharp_flat_plate_mach_05pt00_re1_11pt40e6_isothermal_tw_300_k/README.md) | Sharp flat plate | 5.0 | 11.4e6 | Isothermal 300 K | Falkner-Skan | CFD++ laminar |

## Running a Case

Each case has the same structure:

```
<case_name>/
  data/        reference baseflow data
  figures/     pre-generated comparison figures
  scripts/     validate.py - reproduces all figures
  README.md    case description, conditions, expected outcome
```

To run any case, ensure `simbl` is installed and available, then:

```bash
cd vnv/validation/<solver_name>/<case_name>
python scripts/validate.py
```
