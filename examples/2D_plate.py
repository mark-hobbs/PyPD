"""
Example: 2D plate with a hole
-----------------------------

Run the following command from the root folder:

python -m examples.2D_plate

"""
import numpy as np

import pypd


def build_particle_coordinates(dx, n_div_x, n_div_y):
    """
    Build particle coordinates

    Parameters
    ----------

    Returns
    -------

    Notes
    -----

    """
    particle_coordinates = np.zeros([n_div_x * n_div_y, 2])
    counter = 0

    for i_y in range(n_div_y):  # Depth
        for i_x in range(n_div_x):  # Length
            coord_x = dx * i_x
            coord_y = dx * i_y
            particle_coordinates[counter, 0] = coord_x
            particle_coordinates[counter, 1] = coord_y
            counter += 1

    return particle_coordinates


def build_boundary_conditions(particles, dx):
    bc_flag = np.zeros((len(particles), 2), dtype=np.intc)
    bc_unit_vector = np.zeros((len(particles), 2), dtype=np.intc)

    tol = 1e-6

    for i, particle in enumerate(particles):
        if particle[0] < (0.02 + tol):
            bc_flag[i, 0] = 1
            bc_unit_vector[i, 0] = -1
        if particle[0] > (0.48 - dx - tol):
            bc_flag[i, 0] = 1
            bc_unit_vector[i, 0] = 1

    return bc_flag, bc_unit_vector


def build_hole(particles, centre, radius):
    """
    Build hole in 2D plate

    Parameters
    ----------

    Returns
    -------

    Notes
    -----
    """

    counter = 0
    mask = []

    for particle in particles:
        distance = np.sqrt(
            (particle[0] - centre[0]) ** 2 + (particle[1] - centre[1]) ** 2
        )

        if distance < radius:
            mask.append(counter)

        counter += 1

    return np.delete(particles, mask, axis=0)


def main():
    dx = 1.25e-3
    n_div_x = np.rint(0.5 / dx).astype(int)
    n_div_y = np.rint(0.25 / dx).astype(int)
    hole_centre_x = 0.25 - dx / 2
    hole_centre_y = 0.125 - dx / 2

    x = build_particle_coordinates(dx, n_div_x, n_div_y)
    x = build_hole(x, [hole_centre_x, hole_centre_y], 0.05)
    flag, unit_vector = build_boundary_conditions(x, dx)

    material = pypd.Material(name="quasi-brittle", E=33e9, Gf=130, density=2400, ft=2.5)
    integrator = pypd.EulerCromer()
    bc = pypd.BoundaryConditions(flag, unit_vector, magnitude=1e-4)
    particles = pypd.ParticleSet(x, dx, bc, material)
    linear = pypd.Linear(material, particles, t=dx)
    bonds = pypd.BondSet(particles, linear)
    simulation = pypd.Simulation(dt=None, n_time_steps=20000, damping=0)
    model = pypd.Model(
        particles,
        bonds,
        simulation,
        integrator,
        linear,
    )

    model.run_simulation()
    model.plot_damage(sz=0.75, fig_title="plate-with-a-hole")


main()
