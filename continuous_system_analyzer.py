
from physics_functions import *


# We reserve the symbols 'coor_x', 'coor_y', 'coor_z', as we will be using these for the integrals (only used in the continuous case)

#coor_x, coor_y, coor_z = sp.symbols('coor_x coor_y coor_z')

"""def I_tensor(subscript, coord_s, density_function, boundaries):
    # Subscript is a two-char string ('xx', 'yz' etc...)
    if subscript == 'xx':
        inertia_integrand = function_coordinate_transformation(sp.parse_expr('y**2+z**2'), 'c', coord_s)
    if subscript == 'xy':
        inertia_integrand = function_coordinate_transformation(sp.parse_expr('-x*y'), 'c', coord_s)
    if subscript == 'xz':
        inertia_integrand = function_coordinate_transformation(sp.parse_expr('-x*z'), 'c', coord_s)
    if subscript == 'yx':
        return(I_tensor('xy', coord_s, density_function, boundaries))
    if subscript == 'yy':
        inertia_integrand = function_coordinate_transformation(sp.parse_expr('x**2+z**2'), 'c', coord_s)
    if subscript == 'yz':
        inertia_integrand = function_coordinate_transformation(sp.parse_expr('-y*z'), 'c', coord_s)
    if subscript == 'zx':
        return(I_tensor('xz', coord_s, density_function, boundaries))
    if subscript == 'zy':
        return(I_tensor('yz', coord_s, density_function, boundaries))
    if subscript == 'zz':
        inertia_integrand = function_coordinate_transformation(sp.parse_expr('x**2+y**2'), 'c', coord_s)
    I_res = multidimensional_integral(density_function * inertia_integrand * jacobian[coord_s], boundaries)
    return(sp.simplify(I_res))"""


# -------------------------------------------------------------------------------------
# ----------------------- PHYSICAL CONFIGURATION INPUT --------------------------------
# -------------------------------------------------------------------------------------

# The user selects the coordinate system. The choice will be stored in a variable coord_s in ['c', 'p', 's']

coord_s = ''

print("Input the coordinate system of your choice (empty input defaults to cartesian).")
for key, item in coord_s_descriptions.items():
    print('  ' + item)
while(True):
    coord_s_raw = input(f"  Coordinate system of input displacement vectors: ")
    if coord_s_raw == '':
        coord_s = 'c'
        print("  Cartesian coordinate system selected.")
        break
    elif coord_s_raw.lower() in ['c', 'cart', 'cartesian']:
        coord_s = 'c'
        print("  Cartesian coordinate system selected.")
        break
    elif coord_s_raw.lower() in ['p', 'cyl', 'cylindrical']:
        coord_s = 'p'
        print("  Cylindrical coordinate system selected.")
        break
    elif coord_s_raw.lower() in ['s', 'sph', 'spherical']:
        coord_s = 's'
        print("  Spherical coordinate system selected.")
        break
    print("  Input one of the key strings mentioned in the list of available coordinate systems (e.g. 'c' for cartesian coordinates). Make sure to omit the quotation marks.")


# print out reserved coordinate parameter symbols
coord_symbols = get_reserved_coordinate_parameters(coord_s)
print(f"Reserved symbols in this coordinate system are: {coord_symbols[0]}, {coord_symbols[1]}, {coord_symbols[2]}")

# input three interval boundary conditions
print("Input the boundaries of the physical object as three intervals, using the three reserved symbols and any other free parameters.")
print("Note that the intervals must form a well-ordered dependency set. Input the intervals as two expressions separated by a comma.")
boundary_intervals = False #every item is a list [coordinate symbol, lower bound expression, upper bound expression]. This list will be ordered based on the dependency hierarchy.
while(True):
    
    boundary_intervals = []
    
    for coord_symbol in coord_symbols:
        default_coord_interval = default_coord_intervals[str(coord_symbol)]
        while(True):
            expr_raw = input(f"  Boundary interval of {coord_symbol} (empty=({default_coord_interval[0]},{default_coord_interval[1]})): (")
            if expr_raw == '':
                boundary_intervals.append([coord_symbol, default_coord_interval[0], default_coord_interval[1]])
                break
            try:
                expr_list = expr_raw.split(',')
                boundary_intervals.append([coord_symbol, sp.parse_expr(expr_list[0]), sp.parse_expr(expr_list[1])])
                break
            except IndexError:
                print("  Make sure you input two expressions for each interval - one for each boundary.")
    boundary_intervals = order_boundary_intervals(boundary_intervals)
    if boundary_intervals == False:
        print("  The inputted intervals cannot be ordered in a way that resolves dependency. Try again.")
    else:
        break

# Input the density function
density_func = ""


while(True):
    expr_raw = input(f"  Density scalar field (empty=const. D): ")
    if expr_raw == '':
        density_func = sp.parse_expr("DENSITY")
        break
    try:
        density_func = sp.parse_expr(expr_raw)
        break
    except ValueError:
        print("  Make sure your expression is sanitized.")

#print(density_func)
#print(function_coordinate_transformation(density_func, coord_s, 'c'))


# We now calculate the mass as a triple integral
total_mass = multidimensional_integral(density_func * jacobian[coord_s], boundary_intervals)
print(f"Total mass = {total_mass}")
density_func = density_func.subs(total_mass, "M")
print(f"Density scalar field expressed with the total mass M = {density_func}")


"""

# Identify free symbols and ask user to verify coordinate-dependent assumptions
print("Identifying free symbols and asserting coordinate-based assumptions...")

implicit_free_symbols = []
coordinate_free_symbols = [[], [], []] # matrix of right dependent free symbols; symbols in the last list only appear in the third coordinate, penultimate list symbols only in last two coordinate expressions, first list symbols appear in all three
placed_symbols = []
for i in range(3):
    for object_i in range(number_of_objects):
        for cur_symbol in list(displacements_coord[object_i][i].free_symbols):
            if (not cur_symbol in placed_symbols):
                placed_symbols.append(cur_symbol)
                coordinate_free_symbols[i].append(cur_symbol)


# -------------------------------------------------------------------------------------
# -------------------------------- ASSUMPTIONS ----------------------------------------
# -------------------------------------------------------------------------------------
# We have three kinds of assumptions:
#   1. Automatic assumptions: assumptions made without the user knowing. These currently are: Every parameter is real
#   2. Queried assumptions: assumptions the program asks the user to manually confirm or reads from the user-generated config file
#   3. Conditional assumptions: assumptions that are applied after the final result is printed to inform the user about possible simplifications due to special cases
# Both 2. and 3. are currently in the form of various functions of parameters being assumed positive (like x, 1-x, sin(x) etc...)

# look up which symbols were mentioned in the config file TODO
queried_assumptions_dict = {}
conditional_assumptions_dict = {}

# Assume every symbol is real
assume_real(placed_symbols)
queried_assumptions_dict, conditional_assumptions_dict = get_coord_positive_assumptions(coord_s, coordinate_free_symbols)

print("  Current assumptions:", global_assumptions)
print("  Conditional assumptions:", conditional_assumptions_dict)
"""

#print(function_coordinate_transformation(sp.parse_expr("-x*y"), 'c', coord_s))

# --------------------------------------- COM transformation

com_displacement, com_cartesian_displacement = get_com_displacement_continuous(coord_s, density_func, boundary_intervals, total_mass)

print("Center of mass displacement = ( ", com_displacement[0], ";", com_displacement[1], ";", com_displacement[2], ")")
print("Center of mass displacement cartesian = ( ", com_cartesian_displacement[0], ";", com_cartesian_displacement[1], ";", com_cartesian_displacement[2], ")")

com_density_function = coordinate_shift_continuous(coord_s, density_func, com_displacement)
com_boundary_intervals = coordinate_shift_boundary_intervals(coord_s, boundary_intervals, com_displacement)

print("Density scalar field in the CoM reference frame:")
print(com_density_function)
print("Boundary intervals in the CoM reference frame:")
print(com_boundary_intervals)

# --------------------------------------- MOMENT OF INERTIA TENSOR

# Now we calculate the tensor of inertia for the center of mass

I_com = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
I_com[0][0] = I_tensor_continuous('xx', coord_s, com_density_function, com_boundary_intervals)
I_com[0][1] = I_tensor_continuous('xy', coord_s, com_density_function, com_boundary_intervals)
I_com[0][2] = I_tensor_continuous('xz', coord_s, com_density_function, com_boundary_intervals)
I_com[1][0] = I_tensor_continuous('yx', coord_s, com_density_function, com_boundary_intervals)
I_com[1][1] = I_tensor_continuous('yy', coord_s, com_density_function, com_boundary_intervals)
I_com[1][2] = I_tensor_continuous('yz', coord_s, com_density_function, com_boundary_intervals)
I_com[2][0] = I_tensor_continuous('zx', coord_s, com_density_function, com_boundary_intervals)
I_com[2][1] = I_tensor_continuous('zy', coord_s, com_density_function, com_boundary_intervals)
I_com[2][2] = I_tensor_continuous('zz', coord_s, com_density_function, com_boundary_intervals)

I_com_M = Matrix(I_com.copy())

print("CoM tensor of inertia = (")
#print("  ", I_com[0][0], ";", I_com[0][1], ";", I_com[0][2])
#print("  ", I_com[1][0], ";", I_com[1][1], ";", I_com[1][2])
#print("  ", I_com[2][0], ";", I_com[2][1], ";", I_com[2][2])
print_matrix(I_com)
print(")")




"""

# Now we check whether we need to perform a coordinate transformation

displacements = coordinate_transformation(displacements_coord, coord_s, 'c')
if coord_s != 'c':
    print(f"Displacement vectors after transformation from {coord_s_names[coord_s]} coord. system to cartesian coordinate system:")
    print_vectors(displacements, ('Object n. ', ' displacement = '), 2)

com_displacements = coordinate_shift(displacements, com_displacement)
print("Displacements in the CoM reference frame:")
print_vectors(com_displacements, ('Object n. ', ' displacement = '), 2)
#for i in range(number_of_objects):
#    print(f"  Object n. {i + 1} displacement = ( ", com_displacements[i][0], ";", com_displacements[i][1], ";", com_displacements[i][2], ")")

# Now we calculate the tensor of inertia for the center of mass

I_com = [[0, 0, 0], [0, 0, 0], [0, 0, 0]]
I_com[0][0] = I_tensor('xx', masses, com_displacements)
I_com[0][1] = I_tensor('xy', masses, com_displacements)
I_com[0][2] = I_tensor('xz', masses, com_displacements)
I_com[1][0] = I_tensor('yx', masses, com_displacements)
I_com[1][1] = I_tensor('yy', masses, com_displacements)
I_com[1][2] = I_tensor('yz', masses, com_displacements)
I_com[2][0] = I_tensor('zx', masses, com_displacements)
I_com[2][1] = I_tensor('zy', masses, com_displacements)
I_com[2][2] = I_tensor('zz', masses, com_displacements)

I_com_M = Matrix(I_com.copy())

print("CoM tensor of inertia = (")
#print("  ", I_com[0][0], ";", I_com[0][1], ";", I_com[0][2])
#print("  ", I_com[1][0], ";", I_com[1][1], ";", I_com[1][2])
#print("  ", I_com[2][0], ";", I_com[2][1], ";", I_com[2][2])
print_matrix(I_com)
print(")")

I_com_eigenvectors_raw = I_com_M.eigenvects(simplify=True)
I_com_evec, I_com_eval = unwrap_eigenvectors(I_com_eigenvectors_raw)
PA_stability = [] # Calculates whether rotation about each principal axis is stable

# First check if there's a zero value of I; if yes, that rotation is impossible, and the remaining two are stable
if I_com_eval[0] == 0:
    PA_stability = ['impossible', True, True]
elif I_com_eval[1] == 0:
    PA_stability = [True, 'impossible', True]
elif I_com_eval[2] == 0:
    PA_stability = [True, True, 'impossible']
else:
    for i in range(3):
        PA_stability.append(( (I_com_eval[i] - I_com_eval[(i+1)%3])*(I_com_eval[i] - I_com_eval[(i+2)%3])/(I_com_eval[(i+1)%3]*I_com_eval[(i+2)%3]) >= 0))

def print_I_com_eigenvectors(M_evec, M_eval, V_stability):
    right_strings = []
    for i in range(3):
        if V_stability[i] == 'impossible':
            right_strings.append(f"; I_{i+1} = {M_eval[i]}; Rotation about axis is impossible.")
            #print(f"  e_{i+1} = ({M_evec[i][0]}, {M_evec[i][1]}, {M_evec[i][2]}); I_{i+1} = {M_eval[i]}; Rotation about axis is impossible.")
        else:
            # Try if stability expression is evaluable
            try:
                is_stable = bool(PA_stability[i])
                #print(f"  e_{i+1} = ({M_evec[i][0]}, {M_evec[i][1]}, {M_evec[i][2]}); I_{i+1} = {M_eval[i]}; Rotation about axis is {'stable' if is_stable else 'unstable'}.")
                right_strings.append(f"; I_{i+1} = {M_eval[i]}; Rotation about axis is {'stable' if is_stable else 'unstable'}.")
            except TypeError:
                #print(f"  e_{i+1} = ({M_evec[i][0]}, {M_evec[i][1]}, {M_evec[i][2]}); I_{i+1} = {M_eval[i]}")
                right_strings.append(f"; I_{i+1} = {M_eval[i]}")
    print_vectors(I_com_evec, ('e_', ' = '), 2, 'left', right_strings)


print("Principal axes and their respective eigenvectors:")
print_I_com_eigenvectors(I_com_evec, I_com_eval, PA_stability)
I_com_evec_coord = coordinate_transformation(I_com_evec, 'c', coord_s)
if coord_s != 'c':
    print(f"Principal axes after transformation from cartesian coord. system to {coord_s_names[coord_s]} coordinate system:")
    print_vectors(I_com_evec_coord, ("e'_", ' = '), 2, 'left')

ass_string = full_commit_conditional_assumptions(conditional_assumptions_dict)
print(global_assumptions)
I_com_evec_cond = refine_tensor(I_com_evec_coord)
print(ass_string)
print_vectors(I_com_evec_cond, ("e'_", ' = '), 2, 'left')
"""
