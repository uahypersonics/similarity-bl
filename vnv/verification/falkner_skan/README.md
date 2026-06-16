# Verification — Falkner-Skan solver

Two verification cases that confirm the FS solver recovers known analytic
results. These test the numerical implementation, not the physical model.

## Cases

### 1. Blasius limit (`case_blasius.py`)

At Me = 0.01 (near-incompressible), adiabatic wall, beta = 0 (flat plate),
the compressible Levy-Lees FS system reduces to the classical Blasius equation.

Expected result: f''(0) = 0.469600 (tolerance 1e-3)

Reference: White, Viscous Fluid Flow, 3rd ed., Appendix C.

Note: the Levy-Lees compressible f''(0) at the incompressible limit is 0.4696,
not 0.3321. The two values differ because the Levy-Lees non-dimensionalization
uses a different reference length than the classical Blasius scaling.

### 2. Crocco relation (`case_crocco.py`)

For Pr = 1, adiabatic wall, total enthalpy is conserved across the boundary
layer. This gives the exact pointwise relation (Crocco's theorem):

    tau(eta) = 1 + (gamma - 1)/2 * Me^2 * (1 - fp^2)

Three sub-cases are tested:

| Sub-case | Me  | m   | beta    |
|---|---|---|---|
| A        | 1.5 | 0   | 0.0     |
| B        | 3.0 | 0   | 0.0     |
| C        | 3.0 | 0.1 | 0.1818  |

Sub-case C (non-zero pressure gradient) verifies that the relation holds
independently of beta.

Expected result: max pointwise residual < 1e-3

Reference: Crocco, L. (1932), Sulla trasmissione del calore da una lamina
piana a un fluido scorrente ad alta velocità. Also: White, Viscous Fluid Flow.

## How to run

```bash
cd vnv/verification/falkner_skan
python case_blasius.py
python case_crocco.py
```

Exit code 0 = all PASS, exit code 1 = one or more FAIL.
Crocco figures are saved alongside the script.

## PASS criteria summary

| Case | Quantity | Tolerance |
|---|---|---|
| Blasius | \|f''(0) - 0.4696\| | < 1e-3 |
| Blasius | max\|tau - 1\| | < 1e-3 |
| Crocco A,B,C | max pointwise residual | < 1e-3 |
