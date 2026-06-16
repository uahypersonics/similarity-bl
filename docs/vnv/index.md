# Validation and Verification

These cases confirm that `simbl` produces correct results.

**Validation** compares the similarity solution against independent reference
data (CFD, experiment, published tables) to confirm the model represents
physical reality.

**Verification** checks the numerical implementation against exact or analytic
results to confirm the code is solving the equations correctly.

All scripts and data are in the
[`vnv/`](https://github.com/uahypersonics/similarity-bl/tree/main/vnv)
directory of the repository.

## Verification

*(cases coming soon)*

## Validation

| Case | Geometry | Mach | Re1 (1/m) | Wall BC | Solver | Reference |
|---|---|---|---|---|---|---|
| [Sharp flat plate, M=5, isothermal](validation/sharp_flat_plate_mach_05pt00_re1_11pt40e6_isothermal_tw_300_k/index.md) | Sharp flat plate | 5.0 | 11.4e6 | Isothermal 300 K | Falkner-Skan | CFD++ laminar |
