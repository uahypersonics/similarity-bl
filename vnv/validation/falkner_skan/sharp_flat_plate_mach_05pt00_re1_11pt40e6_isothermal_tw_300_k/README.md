# VnV — Sharp Flat Plate, Mach 5, Isothermal Wall

Validation of the Falkner-Skan similarity solver against a CFD++ laminar baseflow
for a sharp flat plate at Mach 5 with an isothermal wall.

## Flow Conditions

| Parameter | Value |
|---|---|
| Mach number | 5.0 |
| Unit Reynolds number Re₁ | 11.4 × 10⁶ m⁻¹ |
| Freestream temperature T∞ | 75.33 K |
| Wall temperature T_w | 300 K |
| Ratio of specific heats γ | 1.4 |
| Prandtl number Pr | 0.71 |
| Viscosity law | Sutherland |

## Validation Approach

The Falkner-Skan (2D, zero pressure gradient) similarity solution is compared
against wall-normal profiles extracted from the CFD++ laminar baseflow at five
streamwise stations: x = 0.1, 0.2, 0.3, 0.4, 0.5 m.

The similarity profiles are mapped to physical space using the Levy-Lees
transformation for a flat plate (m = 0):

$$\eta_\text{scale} = \sqrt{\frac{\rho_e u_e}{2 \mu_e x}}$$

$$y(\eta) = \frac{1}{\eta_\text{scale}} \int_0^\eta \frac{T}{T_e} \, d\eta'$$

Quantities compared: streamwise velocity ratio u/uₑ, temperature ratio T/Tₑ,
and density ratio ρ/ρₑ.

## Expected Outcome

- **u/uₑ and T/Tₑ**: FS and CFD++ profiles are in close agreement across all
  five stations, confirming that the boundary layer is self-similar and that the
  Levy-Lees transform is applied correctly.
- **ρ/ρₑ**: slight deviation is expected because the CFD++ solution has a weak
  streamwise pressure gradient (ρ/ρₑ ≠ Tₑ/T exactly), while the similarity
  solution assumes constant pressure.

## How to Run

```bash
conda activate dev
cd vnv/sharp_flat_plate_mach_05pt00_re1_11pt40e6_isothermal_tw_300_k
python scripts/validate.py
```

Figures are saved to `figures/`.

## File Layout

```
data/
  base_flow_cfdpp.hdf5          CFD++ baseflow (5 stations, uvel/temp/dens)
figures/
  fs_profiles_eta.png           Similarity profiles in η-space
  similarity_vs_cfdpp_uvel.png  u/uₑ comparison at all x stations
  similarity_vs_cfdpp_temp.png  T/Tₑ comparison at all x stations
  similarity_vs_cfdpp_dens.png  ρ/ρₑ comparison at all x stations
scripts/
  validate.py                   Reproduces all figures
```
