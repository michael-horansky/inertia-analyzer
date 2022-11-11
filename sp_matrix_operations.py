

# Matrix and vector operations which support the sympy library
# These objects are stored as lists (of lists for matrices), where each item is a sympy expression.

# All matrices here represent linear transformations, hence they're all square

# Libraries
import numpy as np
import sympy as sp
from sympy.assumptions import global_assumptions
from sympy.parsing.sympy_parser import parse_expr
from sympy.matrices import Matrix, eye, zeros, ones, diag, GramSchmidt

# -------------------------------------------------------------------------------------
# ----------------------------- TENSOR OPERATIONS -------------------------------------
# -------------------------------------------------------------------------------------

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


# -------------------------------------------------------------------------------------
# ---------------------------- INTEGRAL OPERATIONS ------------------------------------
# -------------------------------------------------------------------------------------

def multidimensional_integral(integrand, boundaries, debug=False):
    # boundaries is a list with element = [param, lower, upper]. They go inwards, so we must flip.
    res = integrand
    for i in range(len(boundaries)-1, -1, -1):
        res = sp.integrate(res, (boundaries[i][0], boundaries[i][1], boundaries[i][2]))
        res = sp.simplify(res)
    if debug:
        debug_str = "  Calculating the integral "
        for i in range(len(boundaries)):
            debug_str += f"S_({boundaries[i][1]})^({boundaries[i][2]}) "
        debug_str += str(integrand)
        for i in range(len(boundaries)-1, -1, -1):
            debug_str += f" d{boundaries[i][0]}"
        debug_str += f" = {res}"
        print(debug_str)
    return(res)



# -------------------------------------------------------------------------------------
# ----------------------- COORDINATE SYSTEM MANAGEMENT --------------------------------
# -------------------------------------------------------------------------------------

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

default_coord_intervals = {
    'x':[-sp.oo, sp.oo],
    'y':[-sp.oo, sp.oo],
    'z':[-sp.oo, sp.oo],
    'rho':[sp.parse_expr("0"), sp.oo],
    'phi':[sp.parse_expr("0"), 2.0 * sp.pi],
    'r':[sp.parse_expr("0"), sp.oo],
    'theta':[sp.parse_expr("0"), sp.pi]
    }

jacobian = {
    'c':sp.parse_expr('1'),
    'p':sp.parse_expr('rho'),
    's':sp.parse_expr('r**2*sin(theta)')
    }

def get_reserved_coordinate_parameters(coord_s):
    res = []
    for param in coord_s_basis[coord_s]:
        res.append(sp.symbols(param))
    return(res)

def coordinate_shift_discrete(vectors_to_shift, shifting_vector):
    # a' = a - v_shift
    result_vectors = []
    for i in range(len(vectors_to_shift)):
        result_vectors.append([vectors_to_shift[i][0] - shifting_vector[0], vectors_to_shift[i][1] - shifting_vector[1], vectors_to_shift[i][2] - shifting_vector[2]])
    return(result_vectors)

def coordinate_shift_continuous(coord_s, func, shifting_vector):
    cartesian_func = function_coordinate_transformation(func, coord_s, 'c')
    cartesian_shifting_vector = tensor_coordinate_transformation(shifting_vector, coord_s, 'c')
    cartesian_func.subs(coord_s_basis[coord_s][0], f"{coord_s_basis[coord_s][0]} - {cartesian_shifting_vector[0]}")
    cartesian_func.subs(coord_s_basis[coord_s][1], f"{coord_s_basis[coord_s][1]} - {cartesian_shifting_vector[1]}")
    cartesian_func.subs(coord_s_basis[coord_s][2], f"{coord_s_basis[coord_s][2]} - {cartesian_shifting_vector[2]}")
    return(function_coordinate_transformation(cartesian_func, 'c', 'coord_s'))

def coordinate_shift_boundary_intervals(coord_s, boundary_intervals, shifting_vector):
    new_boundary_intervals = []
    if coord_s == 'c':
        for i in range(len(boundary_intervals)):
            coord_param_index = coord_s_basis[coord_s].index(str(boundary_intervals[i][0]))
            delta_param = shifting_vector[coord_param_index]
            
            cur_boundary_param, cur_boundary_lower, cur_boundary_upper = boundary_intervals[i][0], boundary_intervals[i][1] - delta_param, boundary_intervals[i][2] - delta_param
            for j in range(3):
                # if the index of the changing coordinate matches the index of the inspected coordinate, we shift the boundaries (which we've already done)
                # otherwise, we substitute in the correct function of the inspected coordinate
                if coord_param_index != j:
                    cur_boundary_lower = cur_boundary_lower.subs(coord_s_basis[coord_s][j], f"{coord_s_basis[coord_s][j]}+{shifting_vector[j]}")
                    cur_boundary_upper = cur_boundary_upper.subs(coord_s_basis[coord_s][j], f"{coord_s_basis[coord_s][j]}+{shifting_vector[j]}")
            new_boundary_intervals.append([cur_boundary_param, cur_boundary_lower, cur_boundary_upper])
    if coord_s == 'p':
        if shifting_vector[0] != 0:
            print("Holy moly! This is a bit too hard!")
            quit()
        for i in range(len(boundary_intervals)):
            # first we shift z' = z - delta z
            # then we sub in z = (z' + delta z) into all the other expressions
            coord_param_index = coord_s_basis[coord_s].index(str(boundary_intervals[i][0]))
            if coord_param_index == 2:
                delta_param = shifting_vector[coord_param_index]
            else:
                delta_param = 0
            new_boundary_intervals.append([boundary_intervals[i][0], boundary_intervals[i][1].subs(coord_s_basis[coord_s][2], f"{coord_s_basis[coord_s][2]}+{shifting_vector[2]}") - delta_param, boundary_intervals[i][2].subs(coord_s_basis[coord_s][2], f"{coord_s_basis[coord_s][2]}+{shifting_vector[2]}") - delta_param])
            #new_boundary_intervals.append([boundary_intervals[i][0], boundary_intervals[i][1] - delta_param, boundary_intervals[i][2] - delta_param])
    if coord_s == 's':
        if shifting_vector[0] != 0:
            print("Holy moly! This is a bit too hard!")
            quit()
        else:
            new_boundary_intervals = boundary_intervals
    return(new_boundary_intervals)

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
    # standardize degenerate cases
    if output_coord_s == 'p' and res[0] == 0:
        res[1] = 0
    if output_coord_s == 's' and res[0] == 0:
        res[1] = 0
        res[2] = 0
    return([sp.refine(sp.simplify(res[0])), sp.refine(sp.simplify(res[1])), sp.refine(sp.simplify(res[2]))])
        
    

def tensor_coordinate_transformation(vectors, input_coord_s, output_coord_s):
    # User-based function. Can pass a vector or a list of vectors.
    # Determine shape by checking if vectors[0] is a list (vector)
    if type(vectors[0]) == list:
        res = []
        for v in vectors:
            res.append(tensor_coordinate_transformation(v, input_coord_s, output_coord_s))
        return(res)
    else:
        return(coordinate_transformation_single_vector(vectors, input_coord_s, output_coord_s))

def coordinate_transformation_single_function(expression, input_coord_s, output_coord_s):
    # here we just sub out reserved coordinate parameters with different ones. EZ!
    #if input_coord_s == output_coord_s:
    #    res = expression
    res = expression
    if input_coord_s == 'c':
        if output_coord_s == 'p':
            # cartesian -> cylindrical
            res = res.subs(coord_s_basis[input_coord_s][0], f"{coord_s_basis[output_coord_s][0]}*cos({coord_s_basis[output_coord_s][1]})")
            res = res.subs(coord_s_basis[input_coord_s][1], f"{coord_s_basis[output_coord_s][0]}*sin({coord_s_basis[output_coord_s][1]})")
        if output_coord_s == 's':
            # cartesian -> spherical
            res = res.subs(coord_s_basis[input_coord_s][0], f"{coord_s_basis[output_coord_s][0]}*sin({coord_s_basis[output_coord_s][1]})*cos({coord_s_basis[output_coord_s][2]})")
            res = res.subs(coord_s_basis[input_coord_s][1], f"{coord_s_basis[output_coord_s][0]}*sin({coord_s_basis[output_coord_s][1]})*sin({coord_s_basis[output_coord_s][2]})")
            res = res.subs(coord_s_basis[input_coord_s][2], f"{coord_s_basis[output_coord_s][0]}*cos({coord_s_basis[output_coord_s][1]})")
            #res = [sp.sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]), sp.atan2(sp.sqrt(v[0]*v[0]+v[1]*v[1]),v[2]), sp.atan2(v[1],v[0])]
    if input_coord_s == 'p':
        if output_coord_s == 'c':
            # cylindrical -> cartesian
            res = res.subs(coord_s_basis[input_coord_s][0], f"sqrt({coord_s_basis[output_coord_s][0]}**2+{coord_s_basis[output_coord_s][1]}**2)")
            res = res.subs(coord_s_basis[input_coord_s][1], f"atan2({coord_s_basis[output_coord_s][1]},{coord_s_basis[output_coord_s][0]})")
            #res = [v[0]*sp.cos(v[1]), v[0]*sp.sin(v[1]), v[2]]
        if output_coord_s == 's':
            # cylindrical -> spherical
            res = res.subs(coord_s_basis[input_coord_s][0], f"{coord_s_basis[output_coord_s][0]}*sin({coord_s_basis[output_coord_s][1]})")
            res = res.subs(coord_s_basis[input_coord_s][2], f"{coord_s_basis[output_coord_s][0]}*cos({coord_s_basis[output_coord_s][1]})")
            #res = [sp.sqrt(v[0]*v[0]+v[2]*v[2]), sp.atan2(v[0],v[2]), v[1]]
    if input_coord_s == 's':
        if output_coord_s == 'c':
            # spherical -> cartesian
            res = res.subs(coord_s_basis[input_coord_s][0], f"sqrt({coord_s_basis[output_coord_s][0]}**2+{coord_s_basis[output_coord_s][1]}**2+{coord_s_basis[output_coord_s][2]}**2)")
            res = res.subs(coord_s_basis[input_coord_s][1], f"atan2(sqrt({coord_s_basis[output_coord_s][0]}**2+{coord_s_basis[output_coord_s][1]}**2),{coord_s_basis[output_coord_s][2]})")
            res = res.subs(coord_s_basis[input_coord_s][2], f"atan2({coord_s_basis[output_coord_s][1]},{coord_s_basis[output_coord_s][0]})")
            
            #res = [v[0] * sp.sin(v[1]) * sp.cos(v[2]), v[0] * sp.sin(v[1]) * sp.sin(v[2]), v[0] * sp.cos(v[1])]
        if output_coord_s == 'p':
            # spherical -> cylindrical
            res = res.subs(coord_s_basis[input_coord_s][0], f"sqrt({coord_s_basis[output_coord_s][0]}**2+{coord_s_basis[output_coord_s][2]}**2)")
            res = res.subs(coord_s_basis[input_coord_s][1], f"atan2({coord_s_basis[output_coord_s][0]},{coord_s_basis[output_coord_s][2]})")
            #res = [v[0] * sp.sin(v[1]), v[2], v[0] * sp.cos(v[1])]
    #if res == expression:
    #    print("Function conversion failed. ({input_coord_s} -> {output_coord_s}) isn't in the supported conversions.")
    #    return(res)
    #return([sp.simplify(res[0]), sp.simplify(res[1]), sp.simplify(res[2])])
    return(sp.simplify(res))

def function_coordinate_transformation(expressions, input_coord_s, output_coord_s):
    # User-based function. Can pass a func or a list of funcs.
    # Determine shape by checking if expressions is a list
    if type(expressions) == list:
        res = []
        for expr in expressions:
            res.append(function_coordinate_transformation(expr, input_coord_s, output_coord_s))
        return(res)
    elif type(expressions) == str:
        return(coordinate_transformation_single_function(sp.parse_expr(expressions), input_coord_s, output_coord_s))
    else:
        return(coordinate_transformation_single_function(expressions, input_coord_s, output_coord_s))

# Boundary intervals ordering function
def order_boundary_intervals(boundary_intervals):
    N = len(boundary_intervals)
    index_symbols = []
    dependency_list = []
    result = []
    for item in boundary_intervals:
        index_symbols.append(item[0])
    for item in boundary_intervals:
        cur_dependency_list = []
        cur_free_symbols = list(item[1].free_symbols) + list(item[2].free_symbols)
        for index_symbol in index_symbols:
            if index_symbol in cur_free_symbols:
                cur_dependency_list.append(index_symbol)
        dependency_list.append([item[0], cur_dependency_list])
    
    # Collapse the dependency list
    while(len(dependency_list)>0):
        # Look for an item with no dependency
        cur_item = False
        for i in range(len(dependency_list)):
            if len(dependency_list[i][1])==0:
                # item's index becomes cur_index, cur_index increases
                result.append(boundary_intervals[index_symbols.index(dependency_list[i][0])])
                cur_item = dependency_list.pop(i)[0]
                break
        # check if we found a dependency-free item
        if cur_item == False:
            return(False)
        # remove dependency
        for i in range(len(dependency_list)):
            if cur_item in dependency_list[i][1]:
                dependency_list[i][1].remove(cur_item)
    return(result)

# -------------------------------------------------------------------------------------
# -------------------------- ASSUMPTIONS MANAGEMENT -----------------------------------
# -------------------------------------------------------------------------------------


def yesno_input(text, default = True, margin = 2):
    print(' '*margin + text, end='')
    default_descriptor = ('yes' if default else 'no')
    while(True):
        raw_answer = input(' [y/n] (empty input=' + default_descriptor + '): ')
        if raw_answer == '':
            return(default)
        if raw_answer.lower() in ['y', 'yes', 'you got it boss', 'google en passant']:
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

# Automatically create assumptions for reserved coordinate parameters
# TODO
def reserved_coord_param_assumptions():
    print("lol")


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





# -------------------------------------------------------------------------------------
# ----------------------------- OUTPUT FUNCTIONS --------------------------------------
# -------------------------------------------------------------------------------------

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


# -------------------------------------------------------------------------------------
# ---------------------------- REFINERY FUNCTIONS -------------------------------------
# -------------------------------------------------------------------------------------
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


# -------------------------------------------------------------------------------------
# --------------------- SYMPY MATRIX OPERATIONS WRAPPERS ------------------------------
# -------------------------------------------------------------------------------------

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




