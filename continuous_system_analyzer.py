
import numpy as np
import sympy as sp
from sympy.parsing.sympy_parser import parse_expr
from sympy.matrices import Matrix, eye, zeros, ones, diag, GramSchmidt

def coordinate_shift(vectors_to_shift, shifting_vector):
    # a' = a - v_shift
    result_vectors = []
    for i in range(len(vectors_to_shift)):
        result_vectors.append([vectors_to_shift[i][0] - shifting_vector[0], vectors_to_shift[i][1] - shifting_vector[1], vectors_to_shift[i][2] - shifting_vector[2]])
    return(result_vectors)

def unwrap_eigenvectors(raw_output, dimension = 3):
    #def general_ordering_func(obj):
    #    return(obj[])
    result = []
    """for i in range(dimension):
        print('Input index', i)
        print(raw_output[i])
        eigenvalue, multiplicity, eigenvector = raw_output[i]
        result.append([eigenvalue, sum(eigenvector[0].tolist(), [])])"""
    for i in range(len(raw_output)):
        eigenvalue, multiplicity, eigenspace = raw_output[i]
        for eigenvector in eigenspace:
            result.append([eigenvalue, sum(eigenvector.tolist(), [])])
    # try if evaluable and sortable
    try:
        # in the general case, we sort by the value of x, then y, then z, then...
        result = sorted(result, key=lambda obj: tuple(obj[1]), reverse=True)
        # 3 dimensional case, we also want to make sure they're right-oriented; we swap last two objects if necessary.
        if dimension == 3:
            # check if sign((e_1 cross e_2)_x) == sign((e_3)_x). If not, swap 2 and 3
            should_swap = False
            if np.sign(result[0][1][1] * result[1][1][2] - result[0][1][2] * result[1][1][1]) != np.sign(result[2][1][0]):
                should_swap = True
            if np.sign(result[0][1][2] * result[1][1][0] - result[0][1][0] * result[1][1][2]) != np.sign(result[2][1][1]):
                should_swap = True
            if np.sign(result[0][1][0] * result[1][1][1] - result[0][1][1] * result[1][1][0]) != np.sign(result[2][1][2]):
                should_swap = True
            if should_swap:
                result[1], result[2] = result[2], result[1]
        return(result)
    except TypeError:
        return(result)

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

# Then the user inputs the discrete system, as a list of mass-displacement tuples which use the declared symbols.

print("Input all point-like masses in the system, in the form [mass,r_x,r_y,r_z]. After the last one, input an empty string.")
masses = []
displacements = [] #d_v[object index][coordinate index]
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
    displacements.append([cur_r_x, cur_r_y, cur_r_z])
    number_of_objects += 1

total_mass = 0
for mass in masses:
    total_mass += mass
total_mass = sp.simplify(total_mass)
print("Total mass = ", total_mass)

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
for i in range(number_of_objects):
    print(f"  Object n. {i + 1} displacement = ( ", com_displacements[i][0], ";", com_displacements[i][1], ";", com_displacements[i][2], ")")

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
print("  ", I_com[0][0], ";", I_com[0][1], ";", I_com[0][2])
print("  ", I_com[1][0], ";", I_com[1][1], ";", I_com[1][2])
print("  ", I_com[2][0], ";", I_com[2][1], ";", I_com[2][2])
print(")")

I_com_eigenvectors_raw = I_com_M.eigenvects(simplify=True)
I_com_ev = unwrap_eigenvectors(I_com_eigenvectors_raw) # I_com_ev[index 0-2] = [eigenvalue, [components, of, eigenvector]]
PA_stability = [] # Calculates whether rotation about each principal axis is stable

# First check if there's a zero value of I; if yes, that rotation is impossible, and the remaining two are stable
if I_com_ev[0][0] == 0:
    PA_stability = ['impossible', True, True]
elif I_com_ev[1][0] == 0:
    PA_stability = [True, 'impossible', True]
elif I_com_ev[2][0] == 0:
    PA_stability = [True, True, 'impossible']
else:
    for i in range(3):
        PA_stability.append(( (I_com_ev[i][0] - I_com_ev[(i+1)%3][0])*(I_com_ev[i][0] - I_com_ev[(i+2)%3][0])/(I_com_ev[(i+1)%3][0]*I_com_ev[(i+2)%3][0]) >= 0))

print("Principal axes and their respective eigenvectors:")
for i in range(3):
    if PA_stability[i] == 'impossible':
        print(f"  e_{i+1} = ({I_com_ev[i][1][0]}, {I_com_ev[i][1][1]}, {I_com_ev[i][1][2]}); I_{i+1} = {I_com_ev[i][0]}; Rotation about axis is impossible.")
    else:
        # Try if stability expression is evaluable
        try:
            is_stable = bool(PA_stability[i])
            print(f"  e_{i+1} = ({I_com_ev[i][1][0]}, {I_com_ev[i][1][1]}, {I_com_ev[i][1][2]}); I_{i+1} = {I_com_ev[i][0]}; Rotation about axis is {'stable' if is_stable else 'unstable'}.")
        except TypeError:
            print(f"  e_{i+1} = ({I_com_ev[i][1][0]}, {I_com_ev[i][1][1]}, {I_com_ev[i][1][2]}); I_{i+1} = {I_com_ev[i][0]}")
#print(sp.simplify(parse_expr(expr_raw)))
