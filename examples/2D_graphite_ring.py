"""
Example: 2D graphite ring under compression

Zhang, X, et al. "Measurement of tensile strength of nuclear graphite based on 
ring compression test." Journal of Nuclear Materials 511 (2018): 134-140.
-------------------------------------------

Run the following command from the root folder:

python -m examples.2D_graphite_ring

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

    return translate_particles(particle_coordinates)


def build_boundary_conditions(particles, dx):
    bc_flag = np.zeros((len(particles), 2), dtype=np.intc)
    bc_unit_vector = np.zeros((len(particles), 2), dtype=np.intc)

    tol = 1e-6

    for i, particle in enumerate(particles):
        if particle[1] < (-0.024 + tol):
            bc_flag[i, 1] = 1
            bc_unit_vector[i, 1] = 1
        if particle[1] > (0.024 - dx - tol):
            bc_flag[i, 1] = 1
            bc_unit_vector[i, 1] = -1

    return bc_flag, bc_unit_vector


def mask_particles_circle(particles, centre, radius, opt):
    """
    Mask particles inside or outside of a circle

    Parameters
    ----------
    opt : str
        Select if particles inside or outside of a circle are deleted
        ["inside" / "outside"]

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

        if opt == "inside":
            if distance <= radius:
                mask.append(counter)
        elif opt == "outside":
            if distance >= radius:
                mask.append(counter)

        counter += 1

    return np.delete(particles, mask, axis=0)


def translate_particles(particles, origin=np.array([0, 0])):
    """
    Translate a set of points to a new origin

    Parameters
    ----------
    origin : str
        Translate a set of points to the origin (a, b)

    Returns
    -------

    Notes
    -----
    https://math.stackexchange.com/questions/1801867/finding-the-centre-of-an-abritary-set-of-points-in-two-dimensions

    """
    centroid = np.average(particles, axis=0)
    return origin - particles + centroid


def main():
    dx = 0.1875e-3  # 3 / 1.5 / 0.75 / 0.375 / 0.1875
    n_div_x = np.rint(0.05 / dx).astype(int)
    n_div_y = np.rint(0.05 / dx).astype(int)
    hole_centre_x = 0.0  # 0.025 - dx/2
    hole_centre_y = 0.0  # 0.025 - dx/2

    x = build_particle_coordinates(dx, n_div_x, n_div_y)
    x = mask_particles_circle(x, [hole_centre_x, hole_centre_y], 0.025 / 2, "inside")
    x = mask_particles_circle(x, [hole_centre_x, hole_centre_y], 0.05 / 2, "outside")
    flag, unit_vector = build_boundary_conditions(x, dx)

    material = pypd.Material(
        name="graphite", E=10e9, Gf=100, density=1780, ft=27.6  # Gf=190
    )
    integrator = pypd.EulerCromer()
    bc = pypd.BoundaryConditions(flag, unit_vector, magnitude=1e-4)
    particles = pypd.ParticleSet(x, dx, bc, material)
    linear = pypd.Linear(material, particles, t=dx)
    bonds = pypd.BondSet(particles, linear)
    simulation = pypd.Simulation(dt=None, n_time_steps=50000, damping=0)
    model = pypd.Model(particles, bonds, simulation, integrator, linear)

    model.run_simulation()
    model.plot_damage(sz=0.5, dsf=0, fig_title="graphite-ring")


main()
