import time

import numpy as np
import sklearn.neighbors as neighbors


def function_timer(func):
    def wrapper(*args, **kwargs):
        start = time.time()
        result = func(*args, **kwargs)
        execution_time = time.time() - start
        print(f"{func.__name__} took {execution_time:.2f} seconds to execute")
        return result

    return wrapper


def build_particle_coordinates(dx, n_div_x, n_div_y, n_div_z):
    """
    TODO: delete this function
    """
    particle_coordinates = np.zeros([n_div_x * n_div_y * n_div_z, 3])
    counter = 0

    for i_z in range(n_div_z):  # Width
        for i_y in range(n_div_y):  # Depth
            for i_x in range(n_div_x):  # Length
                coord_x = dx * i_x
                coord_y = dx * i_y
                coord_z = dx * i_z
                particle_coordinates[counter, 0] = coord_x
                particle_coordinates[counter, 1] = coord_y
                particle_coordinates[counter, 2] = coord_z
                counter += 1

    return particle_coordinates


def build_particle_families(x, horizon):
    """
    Build particle families

    Parameters
    ----------
    x : ndarray (float)
        Material point coordinates in the reference configuration

    horizon : float
        Material point horizon (non-local length scale)

    Returns
    -------
    nlist : list of numpy arrays (int)
        TODO: define a new name and description

    Notes
    -----
    TODO: include a discussion of the algorithm

    """
    n_nodes = np.shape(x)[0]

    tree = neighbors.KDTree(x, leaf_size=160)
    neighbour_list = tree.query_radius(x, r=horizon)

    # Remove identity values, as there is no bond between a node and itself
    neighbour_list = [neighbour_list[i][neighbour_list[i] != i] for i in range(n_nodes)]

    n_family_members = [len(neighbour_list[i]) for i in range(n_nodes)]
    n_family_members = np.array(n_family_members, dtype=np.intc)

    nlist = np.ones((n_nodes, n_family_members.max()), dtype=np.intc) * -1

    for i in range(n_nodes):
        nlist[i, : n_family_members[i]] = neighbour_list[i]

    nlist = nlist.astype(np.intc)

    return nlist, n_family_members


def build_bond_list(nlist):
    """
    Build bond list
    """
    bondlist = [
        [i, j] for i, neighbours in enumerate(nlist) for j in neighbours if i < j
    ]
    return np.array(bondlist, dtype=np.intc)


def build_penetrator():
    """Build rigid penetrator"""
    pass


def build_volume_correction_factors():
    pass


def build_bond_length(x, bondlist):
    """Build the bond length array"""
    n_bonds = np.shape(bondlist)[0]
    xi = np.zeros(
        [
            n_bonds,
        ]
    )
    xi_x = np.zeros(
        [
            n_bonds,
        ]
    )
    xi_y = np.zeros(
        [
            n_bonds,
        ]
    )
    xi_z = np.zeros(
        [
            n_bonds,
        ]
    )

    for k_bond in range(n_bonds):
        node_i = bondlist[k_bond, 0] - 1
        node_j = bondlist[k_bond, 1] - 1

        xi_x[k_bond] = x[node_j, 0] - x[node_i, 0]
        xi_y[k_bond] = x[node_j, 1] - x[node_i, 1]
        xi_z[k_bond] = x[node_j, 2] - x[node_i, 2]
        xi[k_bond] = np.sqrt(xi_x[k_bond] ** 2 + xi_y[k_bond] ** 2 + xi_z[k_bond] ** 2)

    return xi, xi_x, xi_y, xi_z


def calculate_bond_stiffness(E, delta):
    c = (12 * E) / (np.pi * delta**4)
    return c


def calculate_stable_time_step(rho, dx, horizon, c):
    """
    Calculate minimum stable time step

    Parameters
    ----------
    rho : float
        Material density (kg/m^3)

    dx : float
        Mesh resolution (only valid for regular meshes)

    horizon : float
        Horizon radius

    c : float
        Bond stiffness

    Returns
    -------
    dt : float
        Minimum stable time step

    Notes
    -----
    The time step is determined using the stability condition derived in [1]

    [1] Silling, S. A., & Askari, E. (2005). A meshfree method based on the
    peridynamic model of solid mechanics. Computers & structures, 83(17-18),
    1526-1535.

    TODO: is this equation correct and is it valid for 2D and 3D problems?
    """
    return np.sqrt((2 * rho * dx) / (np.pi * horizon**2 * dx * c))
