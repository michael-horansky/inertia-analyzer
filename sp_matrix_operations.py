

# Matrix and vector operations which support the sympy library
# These objects are stored as lists (of lists for matrices), where each item is a sympy expression.

# All matrices here represent linear transformations, hence they're all square

# Libraries
import numpy as np
import sympy as sp
from sympy.assumptions import global_assumptions
from sympy.parsing.sympy_parser import parse_expr
from sympy.matrices import Matrix, eye, zeros, ones, diag, GramSchmidt

# Vector operations

def scalar_product(vec, k):
    result = [0]*len(vec)
    for i in range(len(vec)):
        result[i] = vec[i] * k
    return(result)

def inner_product(vec1, vec2):
    res = 0.0
    for i in range(len(vec1)):
        res += vec1[i] * vec2[i]
    return(res)

def cross_product(vec1, vec2):
    result = [0, 0, 0]
    result[0] = vec1[1] * vec2[2] - vec1[2] * vec2[1] # x component
    result[1] = vec1[2] * vec2[0] - vec1[0] * vec2[2] # y component
    result[2] = vec1[0] * vec2[1] - vec1[1] * vec2[0] # z component
    return(result)

def magnitude(vec):
    return(sp.sqrt(inner_product(vec, vec)))

def normalize_single_vector(vec):
    return(scalar_product(vec, 1.0 / magnitude(vec)))

def normalize(vectors):
    if type(vectors[0]) == list:
        res = []
        for v in vectors:
            res.append(normalize_single_vector(v))
        return(res)
    else:
        return(normalize_single_vector(vectors))

def orthogonal_unit_vector(vec1, vec2):
    return(normalize(cross_product(vec1, vec2)))


# Preset matrices generation

def zero_matrix(d_y, d_x=0):
    res = []
    if d_x == 0:
        d_x = d_y
    for i in range(d_y):
        res.append([])
        for j in range(d_x):
            res[i].append(0)
    return(res)

def identity_matrix(d):
    res = []
    for i in range(d):
        res.append([])
        for j in range(d):
            if i == j:
                res[i].append(1.0)
            else:
                res[i].append(0.0)
    return(res)




# Coordinate system constants and operations

coord_s_names = {
    'c':'cartesian',
    'p':'cylindrical',
    's':'spherical'
    }
coord_s_descriptions = {
    'c': "'c'/'cart'/'cartesian':  cartesian coord. system with parameters x, y, z.",
    'p': "'p'/'cyl'/'cylindrical': cylindrical coord. system with parameters rho (axial dist.), phi (azimuth), z.",
    's': "'s'/'sph'/'spherical':   spherical coord. system with parameters r (radial dist.), theta (inclination), phi (azimuth).",
    }
coord_s_basis = {
    'c':['x', 'y', 'z'],
    'p':['rho', 'phi', 'z'],
    's':['r', 'theta', 'phi']
    }

def yesno_input(text, default = True, margin = 2):
    print(' '*margin + text, end='')
    default_descriptor = ('yes' if default else 'no')
    while(True):
        raw_answer = input(' [y/n] (empty input=' + default_descriptor + '): ')
        if raw_answer == '':
            return(default)
        if raw_answer.lower() in ['y', 'yes', 'you got it boss']:
            return(True)
        if raw_answer.lower() in ['n', 'no' , 'get fucked'     ]:
            return(False)
        print(' '*(margin+2)+"Incomprehensible.", end='')

def assume_real(symbols_list):
    for my_symbol in symbols_list:
        global_assumptions.add(sp.Q.real(my_symbol))
        

# Get queried and conditional assumptions. Queried assumptions are immediately committed, conditionals are stored and returned to apply later

def get_coord_positive_assumptions(coord_s, coord_free_symbols, ask_user = True):
    queried_assumptions_dict = {}
    conditional_assumptions_dict = {}
    if coord_s == 'c':
        return(queried_assumptions_dict, conditional_assumptions_dict)
    elif coord_s == 'p':
        # cylindrical: rho is positive, phi ranges from 0 to 2pi, hence positive
        for rho_like_symbol in coord_free_symbols[0]:
            if yesno_input(f"Can we assume {rho_like_symbol} is positive?", True, 2):
                queried_assumptions_dict[str(rho_like_symbol)] = [rho_like_symbol]
                global_assumptions.add(sp.Q.positive(rho_like_symbol))
            else:
                conditional_assumptions_dict[str(rho_like_symbol)] = [rho_like_symbol]
        for phi_like_symbol in coord_free_symbols[1]:
            if yesno_input(f"Can we assume {phi_like_symbol} is positive?", True, 2):
                queried_assumptions_dict[str(phi_like_symbol)] = [phi_like_symbol]
                global_assumptions.add(sp.Q.positive(phi_like_symbol))
                conditional_assumptions_dict[str(phi_like_symbol)] = [sp.sin(phi_like_symbol), sp.cos(phi_like_symbol)]
            else:
                conditional_assumptions_dict[str(phi_like_symbol)] = [phi_like_symbol, sp.sin(phi_like_symbol), sp.cos(phi_like_symbol)]
    elif coord_s == 's':
        # spherical: r is positive, theta ranges from 0 to pi, hence sin(theta) positive; phi ranges from 0 to 2pi, hence positive
        for r_like_symbol in coord_free_symbols[0]:
            if yesno_input(f"Can we assume {r_like_symbol} is positive?", True, 2):
                queried_assumptions_dict[str(r_like_symbol)] = [r_like_symbol]
                global_assumptions.add(sp.Q.positive(r_like_symbol))
            else:
                conditional_assumptions_dict[str(r_like_symbol)] = [r_like_symbol]
        for theta_like_symbol in coord_free_symbols[1]:
            if yesno_input(f"Can we assume {theta_like_symbol} is in the interval <0, pi>?", True, 2):
                queried_assumptions_dict[str(theta_like_symbol)] = [theta_like_symbol, sp.sin(theta_like_symbol)]
                global_assumptions.add(sp.Q.positive(theta_like_symbol))
                global_assumptions.add(sp.Q.positive(sp.sin(theta_like_symbol)))
            else:
                conditional_assumptions_dict[str(theta_like_symbol)] = [theta_like_symbol, sp.sin(theta_like_symbol)]
        for phi_like_symbol in coord_free_symbols[2]:
            if yesno_input(f"Can we assume {phi_like_symbol} is positive?", True, 2):
                queried_assumptions_dict[str(phi_like_symbol)] = [phi_like_symbol]
                global_assumptions.add(sp.Q.positive(phi_like_symbol))
                conditional_assumptions_dict[str(phi_like_symbol)] = [sp.sin(phi_like_symbol), sp.cos(phi_like_symbol)]
            else:
                conditional_assumptions_dict[str(phi_like_symbol)] = [phi_like_symbol, sp.sin(phi_like_symbol), sp.cos(phi_like_symbol)]
    else:
        print("  Error: coordinate system not supported by get_coord_positive_assumptions. Returning an empty dictionary.")
    return(queried_assumptions_dict, conditional_assumptions_dict)

# Commit conditional assumptions and return descriptor string 

def full_commit_conditional_assumptions(conditional_assumptions_dict):
    res_string = "If we assume that "
    for symbol_name, terms in conditional_assumptions_dict.items():
        for term in terms:
            res_string += str(term) + ", "
            global_assumptions.add(sp.Q.positive(term))
    return(res_string[:-2] + " are positive:")


"""def refine_expression_with_assumptions(expression, positive_assumptions_dict):
    # assumes everything is real, checks what's positive
    cur_symbols = list(expression.free_symbols)
    for cur_symbol in cur_symbols:
"""

def coordinate_shift(vectors_to_shift, shifting_vector):
    # a' = a - v_shift
    result_vectors = []
    for i in range(len(vectors_to_shift)):
        result_vectors.append([vectors_to_shift[i][0] - shifting_vector[0], vectors_to_shift[i][1] - shifting_vector[1], vectors_to_shift[i][2] - shifting_vector[2]])
    return(result_vectors)


def coordinate_transformation_single_vector(v, input_coord_s, output_coord_s):
    # Not to be used by the user.
    res = []
    if input_coord_s == output_coord_s:
        res = v
    if input_coord_s == 'c':
        if output_coord_s == 'p':
            # cartesian -> cylindrical
            res = [sp.sqrt(v[0]*v[0]+v[1]*v[1]), sp.atan2(v[1],v[0]), v[2]]
        if output_coord_s == 's':
            # cartesian -> spherical
            res = [sp.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]), sp.atan2(sp.sqrt(v[0]*v[0]+v[1]*v[1]),v[2]), sp.atan2(v[1],v[0])]
    if input_coord_s == 'p':
        if output_coord_s == 'c':
            # cylindrical -> cartesian
            res = [v[0]*sp.cos(v[1]), v[0]*sp.sin(v[1]), v[2]]
        if output_coord_s == 's':
            # cylindrical -> spherical
            res = [sp.sqrt(v[0]*v[0]+v[2]*v[2]), sp.atan2(v[0],v[2]), v[1]]
    if input_coord_s == 's':
        if output_coord_s == 'c':
            # spherical -> cartesian
            res = [v[0] * sp.sin(v[1]) * sp.cos(v[2]), v[0] * sp.sin(v[1]) * sp.sin(v[2]), v[0] * sp.cos(v[1])]
        if output_coord_s == 'p':
            # spherical -> cylindrical
            res = [v[0] * sp.sin(v[1]), v[2], v[0] * sp.cos(v[1])]
    if res == []:
        print("Vector conversion failed. ({input_coord_s}, {output_coord_s}) isn't in the supported conversions.")
        return(v)
    #return([sp.simplify(res[0]), sp.simplify(res[1]), sp.simplify(res[2])])
    return([sp.refine(res[0]), sp.refine(res[1]), sp.refine(res[2])])
        
    

def coordinate_transformation(vectors, input_coord_s, output_coord_s):
    # User-based function. Can pass a vector or a list of vectors.
    # Determine shape by checking if vectors[0] is a list (vector)
    if type(vectors[0]) == list:
        res = []
        for v in vectors:
            res.append(coordinate_transformation_single_vector(v, input_coord_s, output_coord_s))
        return(res)
    else:
        return(coordinate_transformation_single_vector(vectors, input_coord_s, output_coord_s))



# Printing functions

def print_matrix(matrix, margin=2):
    separator = ' | '
    d_y = len(matrix)
    d_x = len(matrix[0])
    max_width_in_column = [0]*d_x
    offset_matrix = zero_matrix(d_y, d_x)
    for i in range(d_y):
        for j in range(d_x):
            cur_str_len = len(str(matrix[i][j]))
            offset_matrix[i][j] = -cur_str_len
            if cur_str_len > max_width_in_column[j]:
                max_width_in_column[j] = cur_str_len
    breakline_len = sum(max_width_in_column) + (d_x-1)*len(separator)
    breakline = '-'*breakline_len
    for i in range(d_y):
        cur_line = ''
        for j in range(d_x-1):
            left_offset = int(np.floor((offset_matrix[i][j]+max_width_in_column[j])/2))
            cur_line += ' '*left_offset + str(matrix[i][j]) + ' '*((offset_matrix[i][j]+max_width_in_column[j]) - left_offset) + separator
        left_offset = int(np.floor((offset_matrix[i][d_x-1]+max_width_in_column[d_x-1])/2))
        cur_line += ' '*left_offset + str(matrix[i][d_x-1]) + ' '*((offset_matrix[i][d_x-1]+max_width_in_column[d_x-1]) - left_offset)
        print(' '*margin + cur_line)
        if i < d_y-1:
            print(' '*margin + breakline)

def print_vectors(vectors, descriptor = ('', ' = '), margin=2, alignment='left', right_strings = ''):
    d_y = len(vectors)
    d_x = len(vectors[0])
    desc_1, desc_2 = descriptor
    separator = ', '
    line_end = ')'
    if type(right_strings) == str:
        right_strings = [right_strings] * d_y
    max_width_in_column = [0]*d_x
    offset_matrix = zero_matrix(d_y, d_x)
    for i in range(d_y):
        for j in range(d_x):
            cur_str_len = len(str(vectors[i][j]))
            offset_matrix[i][j] = -cur_str_len
            if cur_str_len > max_width_in_column[j]:
                max_width_in_column[j] = cur_str_len
    for i in range(d_y):
        cur_line = desc_1 + str(i+1) + desc_2 + '('
        for j in range(d_x-1):
            if alignment == 'center':
                left_offset = int(np.floor((offset_matrix[i][j]+max_width_in_column[j])/2))
                cur_line += ' '*left_offset + str(vectors[i][j]) + ' '*((offset_matrix[i][j]+max_width_in_column[j]) - left_offset) + separator
            elif alignment == 'left':
                cur_line += str(vectors[i][j]) + ' '*(offset_matrix[i][j]+max_width_in_column[j]) + separator
        if alignment == 'center':
            left_offset = int(np.floor((offset_matrix[i][d_x-1]+max_width_in_column[d_x-1])/2))
            cur_line += ' '*left_offset + str(vectors[i][d_x-1]) + ' '*((offset_matrix[i][d_x-1]+max_width_in_column[d_x-1]) - left_offset) + line_end + right_strings[i]
        elif alignment == 'left':
            cur_line += str(vectors[i][d_x-1]) + ' '*(offset_matrix[i][d_x-1]+max_width_in_column[d_x-1]) + line_end + right_strings[i]
        print(' '*margin + cur_line)


# ------------------------- REFINERY OPERATIONS --------------------
def refine_single_expression(expr):
    return(sp.refine(expr))

def refine_tensor(tensor):
    if type(tensor) != list:
        return(refine_single_expression(tensor))
    else:
        res = []
        for term in tensor:
            res.append(refine_tensor(term))
        return(res)


# sympy matrix operations wrappers

def unwrap_eigenvectors(raw_output, dimension = 3, refine_expr = True):
    # returns a tuple (eigenvectors, eigenvalues)
    result = []
    def unwrap_local_output(output_matrix):
        res_eigenvectors = []
        res_eigenvalues = []
        for i in range(len(output_matrix)):
            if refine_expr:
                res_eigenvalues.append(sp.refine(output_matrix[i][0]))
                res_eigenvectors.append([sp.refine(output_matrix[i][1][0]), sp.refine(output_matrix[i][1][1]), sp.refine(output_matrix[i][1][2])])
            else:
                res_eigenvalues.append(output_matrix[i][0])
                res_eigenvectors.append(output_matrix[i][1])
        return((res_eigenvectors, res_eigenvalues))
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
        return(unwrap_local_output(result))
    except TypeError:
        return(unwrap_local_output(result))




