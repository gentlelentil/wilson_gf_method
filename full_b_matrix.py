def length_a_b(coords_A, coords_B):
    """fct. calculates distance between two coordinates"""
    import numpy as np
    #create AB vector
    vectorAB = np.subtract(coords_A, coords_B)
    #length of AB vector
    length = np.sqrt(np.dot(vectorAB, vectorAB))
    #calculate length: sqrt( (x1 - x2)^2 + (y1 - y2)^2 + (z1 - z2)^2 )
    return length

def angle_a_b_c(coords_A, coords_B, coords_C):
    """fct. calculates the angle between three coordinates"""
    import numpy as np
    #create AB vector
    v1 = np.subtract(coords_A, coords_B)
    #create BC vector
    v2 = np.subtract(coords_C, coords_B)
    #normalise vectors
    #normv1 = np.linalg.norm(v1)
    normv1 = np.sqrt(np.dot(v1, v1))
    normv2 = np.sqrt(np.dot(v2, v2))
    AdotB = np.dot(v1, v2)
    #rounding required for linear systems
    normA_normB= round((normv1 * normv2), 8)
    #calculate angle by arccos((v1 dot v2)/(normv1 * normv2))
    angle = np.arccos((AdotB) / (normA_normB))
    return np.degrees(angle)

def dihedral_a_b_c_d(coords_A, coords_B, coords_C, coords_D):
    import numpy as np
    #create AB vector
    v1 = -1.0 * np.subtract(coords_B, coords_A)
    #create BC vector
    v2 = np.subtract(coords_C, coords_B)
    #norm_v2 = v2 / np.sqrt(np.dot(v2, v2))
    #create CD vector
    v3 = np.subtract(coords_D, coords_C)
    #create plane ABC
    planeABC = 	np.cross(v1, v2)
    #norm_planeABC = np.sqrt(np.dot(planeABC, planeABC))
    #plane BCD
    planeBCD = np.cross(v2, v3)
    #norm_planeBCD = np.sqrt(np.dot(planeBCD, planeBCD))
    planeXplane = np.cross(planeABC, planeBCD)
    y = np.dot(planeXplane, v2)*(1.0 / np.linalg.norm(v2))
    x = np.dot(planeABC, planeBCD)
    dihedral = np.arctan2(y, x)
    #dihedral = np.arccos((np.dot(planeABC, planeBCD)) / (norm_planeABC * norm_planeBCD))
    #check that dihedral doesnt return negative value
    #if np.dot(np.cross(planeABC, planeBCD), norm_v2) < 0:
    #    dihedral = -dihedral
    return np.degrees(dihedral)

def e_vector(A, B):
    import numpy as np
    length = length_a_b(A, B)
    #needs to be vector not length
    norm = np.linalg.norm(A, B)
    e12 = length / norm
    return e12

def length_S_vector(two_atom_list):
    import numpy as np
    e12 = e_vector(two_atom_list[0], two_atom_list[1])
    s_vector = np.append([-e12], [e12])
    return s_vector

def angle_S_vector(three_atom_list):
    import numpy as np
    A, B, C = three_atom_list[:2]
    angle = angle_a_b_c(A, B, C)
    length_1_2 = length_a_b(A, B)
    length_2_3 = length_a_b(B, C)
    e_vectors = []
    for id_a, a in enumerate(three_atom_list):
        for id_b, b in enumerate(three_atom_list):
            if id_a == id_b:
                continue
            e_vector = e_vector(a, b)
            if id_a > id_b:
                e_vector = -e_vector
            e_vectors.append(e_vector)
        
    s1 = (np.cos(angle) * (e_vectors[4] - e_vectors[5]) / (length_1_2 * np.sin(angle)))
    s2 = (np.cos(angle) * (e_vectors[5] - e_vectors[4])) / (length_2_3 * np.sin(angle))
    s3 = ((length_1_2 - length_2_3 * np.cos(angle) * e_vectors[4] + (length_2_3 - length_1_2 * np.cos(angle) * e_vectors[5]))) / (length_1_2 * length_2_3 * np.sin(angle))
    s_vector = np.concatenate((s1, s2, s3))
    return s_vector

def dihedral_S_vector(four_atom_list):
    import numpy as np
    A, B, C, D = four_atom_list[:4]
    angle123 = angle_a_b_c(A, B, C)
    angle234 = angle_a_b_c(B, C, D)   
    length_1_2 = length_a_b(A, B)
    length_2_3 = length_a_b(B, C)
    length_3_4 = length_a_b(C, D)
    e_vectors = []
    for id_a, a in enumerate(four_atom_list):
        for id_b, b in enumerate(four_atom_list):
            if id_a == id_b:
                continue
            e_vector = e_vector(a, b)
            if id_a > id_b:
                e_vector = -e_vector
                #consider if necessary, maybe coordinates other way around produce negative e vector naturally
            e_vectors.append(e_vector)
    s1 = (-1 * (np.cross(e_vectors[0], e_vectors[4])) / (length_1_2*(np.sin(angle123)**2)))
    s2 = ((length_2_3 - length_1_2 * np.cos(angle123) / (length_2_3 * length_1_2 * np.sin(angle123))) * ((np.cross(e_vectors[0], e_vectors[4])) / np.sin(angle123)) + (np.cos(angle234) / (length_2_3 * np.sin(angle234))) * ((np.cross(e_vectors[11], e_vectors[7])) / np.sin(angle234)))
    s3 = ((length_2_3 - length_3_4 * np.cos(angle234) / (length_2_3 * length_3_4 * np.sin(angle234))) * ((np.cross(e_vectors[11], e_vectors[7])) / np.sin(angle234)) + (np.cos(angle123) / (length_2_3 * np.sin(angle123))) * ((np.cross(e_vectors[0], e_vectors[4])) / np.sin(angle123)))
    s4 = (-1 * (np.cross(e_vectors[11], e_vectors[7])) / (length_3_4*(np.sin(angle234)**2)))
    s_vector = np.concatenate((s1, s2, s3, s4))
    return s_vector

def full_B_matrix(atom_coord_list):
    import numpy as np
    if len(atom_coord_list) < 2:
        return 0
    num_s_vectors = list(3 * len(atom_coord_list) - 6)
    if len(atom_coord_list) == 2:
        num_s_vectors = 1
    s_vectors = []
    for id, svec in enumerate(num_s_vectors):
        if svec == 1:
            array = np.zeros((1, len(atom_coord_list)))
            twoatomlist = atom_coord_list[:2]
            s_vec = length_S_vector(twoatomlist)
            array[0] = s_vec[0]
            array[1] = s_vec[1]
            s_vectors.append(array)
        if svec == 2:
            array = np.zeros((1, len(atom_coord_list)))
            threeatomlist = atom_coord_list[:3]
            array[0] = s_vec[0]
            array[1] = s_vec[1]
            array[2] = s_vec[2]
            s_vectors.append(array)
        if svec == 3:
            array = np.zeros((1, len(atom_coord_list)))
            twoatomlist = atom_coord_list[0:2]
            s_vec = length_S_vector(twoatomlist)
            array[1] = s_vec[0]
            array[2] = s_vec[1]
            s_vectors.append(array)
        if svec -1 % 3 == 0:
            #angles
            array = np.zeros((1, len(atom_coord_list)))
            initial_index = (id / 3) + 2
            threeatomlist = atom_coord_list[initial_index -2 : initial_index]
            s_vec = angle_S_vector(threeatomlist)
            array[initial_index - 2] = s_vec[0]
            array[initial_index - 1] = s_vec[1]
            array[initial_index] = s_vec[2]
            s_vectors.append(array)
        if svec - 2 % 3 == 0:
            #bonds
            array = np.zeros((1, len(atom_coord_list)))
            initial_index = ((id + 2) / 3) + 1
            twoatomlist = atom_coord_list[initial_index - 2: initial_index - 1]
            s_vec = length_S_vector(twoatomlist)
            array[initial_index - 2] = s_vec[0]
            array[initial_index - 1] = s_vec[1]
            s_vectors.append(array)
        if svec % 3 == 0:
            #dihedrals
            initial_index = (id + 1 / 3) + 1
            fouratomlist = atom_coord_list[initial_index - 3: initial_index]
            s_vec = dihedral_S_vector(fouratomlist)
            array[initial_index - 3] = s_vec[0]
            array[initial_index - 2] = s_vec[1]
            array[initial_index - 1] = s_vec[2]
            array[initial_index] = s_vec[3]
            s_vectors.append(array)
    bmat = np.stack(s_vectors[:])
    return bmat        