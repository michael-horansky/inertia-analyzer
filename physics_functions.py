
from sp_matrix_operations import *



# -------------------- Center of mass ---------------------

def get_com_displacement_discrete(masses, displacements, total_mass):
    com_displacement = [0, 0, 0]
    for i in range(len(masses)):
        com_displacement[0] += masses[i] * displacements[i][0]
        com_displacement[1] += masses[i] * displacements[i][1]
        com_displacement[2] += masses[i] * displacements[i][2]
    com_displacement[0] = sp.simplify(com_displacement[0] / total_mass)
    com_displacement[1] = sp.simplify(com_displacement[1] / total_mass)
    com_displacement[2] = sp.simplify(com_displacement[2] / total_mass)
    return(com_displacement)

def get_com_displacement_continuous(coord_s, density_function, boundaries, total_mass):
    """if coord_s == 'c':
        com_x = multidimensional_integral(density_function * sp.parse_expr('x') * jacobian[coord_s], boundaries) / total_mass
        com_y = multidimensional_integral(density_function * sp.parse_expr('y') * jacobian[coord_s], boundaries) / total_mass
        com_z = multidimensional_integral(density_function * sp.parse_expr('z') * jacobian[coord_s], boundaries) / total_mass
        return([com_x, com_y, com_z])
    if coord_s == 'p':
        com_rho = multidimensional_integral(density_function * sp.parse_expr('rho') * jacobian[coord_s], boundaries) / total_mass
        com_phi = sp.acos(multidimensional_integral(density_function * sp.parse_expr('cos(theta)') * jacobian[coord_s], boundaries)) / total_mass
        com_z   = multidimensional_integral(density_function * sp.parse_expr('z') * jacobian[coord_s], boundaries) / total_mass
        return([com_rho, com_phi, com_z])
    if coord_s == 's':
        print(density_function * sp.parse_expr('r') * jacobian[coord_s])
        com_r     = multidimensional_integral(density_function * sp.parse_expr('r') * jacobian[coord_s], boundaries) / total_mass
        com_theta = sp.asin(multidimensional_integral(density_function * sp.parse_expr('sin(theta)') * jacobian[coord_s], boundaries)) / total_mass
        com_phi   = sp.acos(multidimensional_integral(density_function * sp.parse_expr('cos(phi)') * jacobian[coord_s], boundaries)) / total_mass
        return([com_r, com_theta, com_phi])"""
    integrand_vector = function_coordinate_transformation(['x', 'y', 'z'], 'c', coord_s)
    com_cartesian_vector = [0, 0, 0]
    com_cartesian_vector[0] = multidimensional_integral(density_function * integrand_vector[0] * jacobian[coord_s], boundaries) / total_mass
    com_cartesian_vector[1] = multidimensional_integral(density_function * integrand_vector[1] * jacobian[coord_s], boundaries) / total_mass
    com_cartesian_vector[2] = multidimensional_integral(density_function * integrand_vector[2] * jacobian[coord_s], boundaries) / total_mass
    com_og_vector = tensor_coordinate_transformation(com_cartesian_vector, 'c', coord_s)
    return(com_og_vector, com_cartesian_vector)


# --------------------- Moment of inertia tensor -------------------

def I_tensor_discrete(subscript, masses, disps):
    # Subscript is a two-char string ('xx', 'yz' etc...)
    I_res = 0
    if subscript == 'xx':
        for i in range(len(masses)):
            I_res += masses[i] * (disps[i][1] * disps[i][1] + disps[i][2] * disps[i][2])
    if subscript == 'xy':
        for i in range(len(masses)):
            I_res -= masses[i] * disps[i][0] * disps[i][1]
    if subscript == 'xz':
        for i in range(len(masses)):
            I_res -= masses[i] * disps[i][0] * disps[i][2]
    if subscript == 'yx':
        return(I_tensor_discrete('xy', masses, disps))
    if subscript == 'yy':
        for i in range(len(masses)):
            I_res += masses[i] * (disps[i][0] * disps[i][0] + disps[i][2] * disps[i][2])
    if subscript == 'yz':
        for i in range(len(masses)):
            I_res -= masses[i] * disps[i][1] * disps[i][2]
    if subscript == 'zx':
        return(I_tensor_discrete('xz', masses, disps))
    if subscript == 'zy':
        return(I_tensor_discrete('yz', masses, disps))
    if subscript == 'zz':
        for i in range(len(masses)):
            I_res += masses[i] * (disps[i][0] * disps[i][0] + disps[i][1] * disps[i][1])
    return(sp.simplify(I_res))

def I_tensor_continuous(subscript, coord_s, density_function, boundaries):
    # Subscript is a two-char string ('xx', 'yz' etc...)
    if subscript == 'xx':
        inertia_integrand = function_coordinate_transformation(sp.parse_expr('y**2+z**2'), 'c', coord_s)
    if subscript == 'xy':
        inertia_integrand = function_coordinate_transformation(sp.parse_expr('-x*y'), 'c', coord_s)
    if subscript == 'xz':
        inertia_integrand = function_coordinate_transformation(sp.parse_expr('-x*z'), 'c', coord_s)
    if subscript == 'yx':
        return(I_tensor_continuous('xy', coord_s, density_function, boundaries))
    if subscript == 'yy':
        inertia_integrand = function_coordinate_transformation(sp.parse_expr('x**2+z**2'), 'c', coord_s)
    if subscript == 'yz':
        inertia_integrand = function_coordinate_transformation(sp.parse_expr('-y*z'), 'c', coord_s)
    if subscript == 'zx':
        return(I_tensor_continuous('xz', coord_s, density_function, boundaries))
    if subscript == 'zy':
        return(I_tensor_continuous('yz', coord_s, density_function, boundaries))
    if subscript == 'zz':
        inertia_integrand = function_coordinate_transformation(sp.parse_expr('x**2+y**2'), 'c', coord_s)
    I_res = multidimensional_integral(density_function * inertia_integrand * jacobian[coord_s], boundaries, debug=True)
    return(sp.simplify(I_res))
