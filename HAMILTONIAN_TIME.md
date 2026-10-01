# Hamiltonian time and prescribed equilibrium potential

Hamiltonian-time polynomial coefficients include the electric-drift curl
in the strong-electric-field symplectic Jacobian. The parallel-current
coefficient uses the same Jacobian. The weak-field coefficients retain
only the magnetic terms.

A configured `er_profile_file` now supplies the potential in both weak and
strong electric-field orbit models. This keeps the prescribed potential
fixed when changing the orbit ordering. Weak ordering does not add the
strong-field electric-drift kinetic or symplectic terms. Without an Er
profile, the weak branch retains its legacy `eps_Phi*A_phi` potential.
An explicitly retained externally supplied potential is preserved.

The Hamiltonian-coefficient test independently constructs the local
symplectic form and checks both particle-charge signs and field orderings.
Full density and field response require additional validation in the
response applet.
