
from sp_matrix_operations import *


# We reserve the symbols 'coor_x', 'coor_y', 'coor_z', as we will be using these for the integrals (only used in the continuous case)

coor_x, coor_y, coor_z = sp.symbols('coor_x coor_y coor_z')

def I_tensor(subscript, masses, disps):
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
        return(I_tensor('xy', masses, disps))
    if subscript == 'yy':
        for i in range(len(masses)):
            I_res += masses[i] * (disps[i][0] * disps[i][0] + disps[i][2] * disps[i][2])
    if subscript == 'yz':
        for i in range(len(masses)):
            I_res -= masses[i] * disps[i][1] * disps[i][2]
    if subscript == 'zx':
        return(I_tensor('xz', masses, disps))
    if subscript == 'zy':
        return(I_tensor('yz', masses, disps))
    if subscript == 'zz':
        for i in range(len(masses)):
            I_res += masses[i] * (disps[i][0] * disps[i][0] + disps[i][1] * disps[i][1])
    return(sp.simplify(I_res))


# First the user declares all parameters of the system

"""
symbol_names_raw = input("Input the list of variable parameters, as a string of names separated by commas: ")

symbol_names = symbol_names_raw.split(' ')
parameters = []
for name in symbol_names:
    parameters.append(sp.symbols(name))"""


#print_matrix([['a', 'aaa', 'aaaaa'],['aaaaa', 'aa', 'aaaa'],['aaaaa', 'aaa', 'a']])


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


# converting I_ij between coordinates: we have the basis transformation matrix B_ij', B_ij' w_i' = w_i.
# Then L_ij = I_ij (B_ij' w_i') = (I_ij B_ij') w_i' (by associativity) = I_ij' w_i', I_ij' = I_ij B_ij'

# Then the user inputs the discrete system, as a list of mass-displacement tuples which use the declared symbols.

print(f"Input all point-like masses in the system, in the form [mass,r_{coord_s_basis[coord_s][0]},r_{coord_s_basis[coord_s][1]},r_{coord_s_basis[coord_s][2]}]. After the last one, input an empty string.")
masses = []
displacements = [] #d_v[object index][coordinate index]
displacements_coord = [] #d_v[object index][coordinate index]
number_of_objects = 0

while(True):
    expr_raw = input(f"  Parameters of object n. {number_of_objects + 1}: ")
    if expr_raw == '':
        break
    expr_list = expr_raw.split(',')
    try:
        cur_mass = parse_expr(expr_list[0])
        cur_r_x, cur_r_y, cur_r_z = parse_expr(expr_list[1]), parse_expr(expr_list[2]), parse_expr(expr_list[3])
    except IndexError:
        print("  Make sure you input four expressions for each object.")
        continue
    masses.append(cur_mass)
    displacements_coord.append([cur_r_x, cur_r_y, cur_r_z])
    number_of_objects += 1


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


# ------------------------- ASSUMPTIONS ---------------------------
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


total_mass = 0
for mass in masses:
    total_mass += mass
total_mass = sp.simplify(total_mass)
print("Total mass = ", total_mass)

# Now we check whether we need to perform a coordinate transformation

displacements = coordinate_transformation(displacements_coord, coord_s, 'c')
if coord_s != 'c':
    print(f"Displacement vectors after transformation from {coord_s_names[coord_s]} coord. system to cartesian coordinate system:")
    print_vectors(displacements, ('Object n. ', ' displacement = '), 2)

# Now we find the center of mass, so we can move there.

com_displacement = [0, 0, 0]
for i in range(number_of_objects):
    com_displacement[0] += masses[i] * displacements[i][0]
    com_displacement[1] += masses[i] * displacements[i][1]
    com_displacement[2] += masses[i] * displacements[i][2]
com_displacement[0] = sp.simplify(com_displacement[0] / total_mass)
com_displacement[1] = sp.simplify(com_displacement[1] / total_mass)
com_displacement[2] = sp.simplify(com_displacement[2] / total_mass)

print("Center of mass displacement = ( ", com_displacement[0], ";", com_displacement[1], ";", com_displacement[2], ")")

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

"""for i in range(3):
    if PA_stability[i] == 'impossible':
        print(f"  e_{i+1} = ({M_evec[i][0]}, {M_evec[i][1]}, {M_evec[i][2]}); I_{i+1} = {M_eval[i]}; Rotation about axis is impossible.")
    else:
        # Try if stability expression is evaluable
        try:
            is_stable = bool(PA_stability[i])
            print(f"  e_{i+1} = ({M_evec[i][0]}, {M_evec[i][1]}, {M_evec[i][2]}); I_{i+1} = {M_eval[i]}; Rotation about axis is {'stable' if is_stable else 'unstable'}.")
        except TypeError:
            print(f"  e_{i+1} = ({M_evec[i][0]}, {M_evec[i][1]}, {M_evec[i][2]}); I_{i+1} = {M_eval[i]}")"""

#print(sp.simplify(parse_expr(expr_raw)))
