# Molecular dynamics simulation (LAMMPS) @ Coulomb crystal $^{40}Ca^+$.

# Ca⁺ Single-Species Cooling (PyLiON + LAMMPS)

This repository contains a Python script that simulates **Ca⁺ ion cooling** in a linear Paul trap using **PyLiON** as a front-end to **LAMMPS** (`units si`). The script:

1. Builds a compact cloud of Ca⁺ ions.
2. Configures a linear Paul trap from a target Mathieu parameter \( q \).
3. Pre-minimizes in the pseudopotential (no micromotion).
4. Thermalizes with a Langevin bath at 1 mK.
5. Monitors the system with time-averaged velocities and dumps.
6. Applies staged laser cooling to Ca⁺ and then coasts with no baths.
7. Patches the LAMMPS input file (boundaries, timestep, etc.).
8. Runs LAMMPS, analyzes the final structure, and generates plots/CSV logs.

All user-level inputs are in **SI units** (m, s, kg, C, V), and the final plots are in **µm** for readability.

---

## 1. Physical Model

We simulate \( N = 70 \) identical Ca⁺ ions:

- Charge: \( Q = +e \)
- Mass: \( m \approx 40 \,\text{amu} \)

Each ion \( i \) has position \( \mathbf{r}_i = (x_i, y_i, z_i) \) and velocity \( \mathbf{v}_i \). The equations of motion are

\[
m \frac{d^2 \mathbf{r}_i}{dt^2}
= \mathbf{F}_{\text{trap}}(\mathbf{r}_i, t)
+ \mathbf{F}_{\text{Coul}}^{(i)}
+ \mathbf{F}_{\text{Langevin}}^{(i)}
+ \mathbf{F}_{\text{laser}}^{(i)}.
\]

LAMMPS integrates these using a velocity–Verlet scheme with a timestep

\[
\Delta t = 5\times 10^{-10} \,\text{s}.
\]

---

## 2. Linear Paul Trap and Mathieu Parameter

The ions are confined in a **linear Paul trap**. An idealized quadrupole potential is

\[
\Phi(x,y,z,t) =
\frac{V_{\text{RF}}\cos(\Omega t)}{2 r_0^2}(x^2 - y^2)
+ \Phi_{\text{DC}}(z),
\]

where:

- \( V_{\text{RF}} \): RF voltage amplitude  
- \( \Omega = 2\pi f_{\text{RF}} \): RF angular frequency  
- \( r_0 \): characteristic electrode radius  
- \( \Phi_{\text{DC}}(z) \): axial DC potential

In the script:

- \( r_0 = 1\, \text{mm} \)
- \( f_{\text{RF}} = 5.634 \times 10^6\, \text{Hz} \) (5.634 MHz)
- Axial DC confinement from endcaps: `endcapvoltage = 20.0` V
- Geometry factor: \( \kappa = 0.3 \)
- Axial length scale: `length = 3.5e-3` m

The radial motion in a Paul trap is described by a **Mathieu equation** with radial parameter

\[
q_r = \frac{2 Q V_{\text{RF}}}{m r_0^2 \Omega^2}.
\]

We choose a target value \( q_{\text{target}} = 0.25 \) and solve for \( V_{\text{RF}} \) using PyLiON’s

```python
Vrf, _ = pl.trapaqtovoltage(ca_ions, trap, 0, q_target)
trap['voltage'] = abs(Vrf)

