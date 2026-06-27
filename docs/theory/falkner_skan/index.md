# Falkner-Skan Equations

The Falkner-Skan (FS) equations are a obtained by reducing the [2D compressible boundary layer equations](../boundary_layer_equations/2d_equations.md) to an ODE system. All assumptions from the boundary layer equations carry over.

??? abstract "Inherited assumptions"
    - Steady, two-dimensional, laminar flow
    - Calorically perfect gas ($\gamma$, $c_p$, $\mathrm{Pr}$ constant)
    - Temperature-dependent viscosity $\mu = \mu(T)$
    - Thin boundary layer ($\delta \ll L$)

## Similarity Ansatz

$$u_e = C x^m, \qquad \psi = \sqrt{2\xi}\,f(\eta), \qquad u = u_e\,f'(\eta)$$

## Outer Flow

$$-\frac{dp}{dx} = \rho_e u_e \frac{du_e}{dx}$$

## Definitions

**Stream function** [^white_2006][^schlichting_2017]:

$$\rho u = \frac{\partial\psi}{\partial y}, \qquad \rho v = -\frac{\partial\psi}{\partial x}$$

**Levy-Lees similarity coordinates** [^choen_1955][^white_2006]:

$$\xi = \int_0^x \rho_e \mu_e u_e\,dx', \qquad \eta = \frac{u_e}{\sqrt{2\xi}}\int_0^y \rho\,dy'$$

**Hartree parameter**, **Chapman-Rubesin factor**, **temperature ratio**:

$$\beta_H = \frac{2m}{m+1}, \qquad C = \frac{\rho\mu}{\rho_e\mu_e}, \qquad \tau = \frac{T}{T_e}$$

**Edge Mach number**, **Prandtl number**:

$$M_e = \frac{u_e}{\sqrt{\gamma R T_e}}, \qquad \mathrm{Pr} = \frac{\mu c_p}{k}$$

## ODE System

The [2D compressible BL equations](../boundary_layer_equations/2d_equations.md) reduce to ODEs in $\eta$ (see derivation below).

**x-momentum**:

$$(Cf'')' + ff'' + \beta_H(\tau - f'^2) = 0$$

**Energy**:

$$\left(\frac{C}{\mathrm{Pr}}\tau'\right)' + f\tau' + (\gamma-1)M_e^2\!\left[Cf''^2 - \beta_H\tau f'\right] = 0$$

## Boundary Conditions

**Wall** ($\eta = 0$):

$$f(0) = 0, \qquad f'(0) = 0$$

- **Isothermal**: $\tau(0) = T_w/T_e$ (prescribed)
- **Adiabatic**: $\tau'(0) = 0$

**Edge** ($\eta \to \infty$):

$$f' = 1, \qquad \tau = 1$$

[^choen_1955]: Cohen, C. B. & Reshotko, E. (1955). *Similar solutions for the compressible laminar boundary layer with heat transfer and pressure gradient*. NACA TN 1293. [PDF](https://apps.dtic.mil/sti/tr/pdf/ADA379809.pdf)
[^white_2006]: White, F. M. (2006). *Viscous Fluid Flow*, 3rd ed. McGraw-Hill, New York.
[^schlichting_2017]: Schlichting, H. & Gersten, K. (2017). *Boundary Layer Theory*, 9th ed. Springer. DOI: [10.1007/978-3-662-52919-5](https://doi.org/10.1007/978-3-662-52919-5)

