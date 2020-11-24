
def length_a_b(coords_A, coords_B):
    """fct. calculates distance between two coordinates"""
    import numpy as np
    vectorAB = np.subtract(coords_A, coords_B)
    length = np.sqrt(np.dot(vectorAB, vectorAB))
    return length

def angle_a_b_c(coords_A, coords_B, coords_C):
    """fct. calculates the angle between three coordinates"""
    import numpy as np
    v1 = coords_A - coords_B
    v2 = coords_C - coords_B
    lengthv1v2 = length_a_b(v1, v2)
    angle = np.arccos(np.dot(v1, v2) / lengthv1v2)
    return angle

def e_vector(A, B):
    import numpy as np
    vector = B - A
    length = length_a_b(A, B)
    e12 = vector / length
    return e12

def length_S_vector(two_atom_list):
    import numpy as np
    e12 = e_vector(two_atom_list[0], two_atom_list[1])
    s_vector = np.append([-e12], [e12])
    return s_vector

def angle_S_vector(three_atom_list):
    import numpy as np
    A, B, C = three_atom_list[:3]
    angle = angle_a_b_c(A, B, C)
    length_1_2 = length_a_b(A, B)
    length_2_3 = length_a_b(B, C)
    e_vectors = []
    for id_a, a in enumerate(three_atom_list):
        for id_b, b in enumerate(three_atom_list):
            if id_a == id_b:
                continue
            e_vec = e_vector(a, b)
            e_vectors.append(e_vec)
    s1 = ((np.cos(angle) * e_vectors[2]) - e_vectors[3]) / (length_1_2 * np.sin(angle))
    s3 = ((np.cos(angle) * e_vectors[3]) - e_vectors[2]) / (length_2_3 * np.sin(angle))
    #s2 = ((length_1_2 - length_2_3 * np.cos(angle) * e_vectors[2] + (length_2_3 - length_1_2 * np.cos(angle) * e_vectors[3]))) / (length_1_2 * length_2_3 * np.sin(angle))
    s2 = -1 * (s1 + s3)
    s_vector = np.concatenate((s1, s2, s3), axis=None)
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
            e_vec = e_vector(a, b)
            e_vectors.append(e_vec)
    s1 = (-1 * (np.cross(e_vectors[0], e_vectors[4])) / (length_1_2*(np.sin(angle123)**2)))
    s2 = ((length_2_3 - length_1_2 * np.cos(angle123) / (length_2_3 * length_1_2 * np.sin(angle123))) * ((np.cross(e_vectors[0], e_vectors[4])) / np.sin(angle123)) + (np.cos(angle234) / (length_2_3 * np.sin(angle234))) * ((np.cross(e_vectors[11], e_vectors[7])) / np.sin(angle234)))
    s3 = ((length_2_3 - length_3_4 * np.cos(angle234) / (length_2_3 * length_3_4 * np.sin(angle234))) * ((np.cross(e_vectors[11], e_vectors[7])) / np.sin(angle234)) + (np.cos(angle123) / (length_2_3 * np.sin(angle123))) * ((np.cross(e_vectors[0], e_vectors[4])) / np.sin(angle123)))
    s4 = (-1 * (np.cross(e_vectors[11], e_vectors[7])) / (length_3_4*(np.sin(angle234)**2)))
    s_vector = np.concatenate((s1, s2, s3, s4), axis=None)
    return s_vector

def full_B_matrix(atom_coord_list):
    import numpy as np
    if len(atom_coord_list) < 2:
        print("No reflections possible for one atom. At least two are required")
        return 0
    tot_s_vectors = 3 * len(atom_coord_list) - 6
    print("Number of atoms found:", len(atom_coord_list))
    print("Minimum DOF for", len(atom_coord_list), "atoms are:", tot_s_vectors)
    num_s_vectors = range(1, tot_s_vectors + 1)
    if len(atom_coord_list) == 2:
        num_s_vectors = [1]
    s_vectors = []
    for id, svec in enumerate(num_s_vectors):
        if svec == 1:
            array = np.zeros((len(atom_coord_list) * 3))
            twoatomlist = atom_coord_list[:2]
            s_vec = length_S_vector(twoatomlist)
            array[0:3] = s_vec[0:3] 
            array[3:6] = s_vec[3:6]
            s_vectors.append(array)
            continue
        if svec == 2:
            array = np.zeros((len(atom_coord_list) * 3))
            threeatomlist = atom_coord_list[:3]
            s_vec = angle_S_vector(threeatomlist)
            array[0:3] = s_vec[0:3]
            array[3:6] = s_vec[3:6]
            array[6:9] = s_vec[6:]
            s_vectors.append(array)
            continue
        if svec == 3:
            array = np.zeros((len(atom_coord_list) * 3))
            twoatomlist = atom_coord_list[0:2]
            s_vec = length_S_vector(twoatomlist)
            array[3:6] = s_vec[0:3]
            array[6:9] = s_vec[3:]
            s_vectors.append(array)
            continue
        if (svec - 1) % 3 == 0:
            #angles
            array = np.zeros((len(atom_coord_list) * 3))
            atom_index = int(((svec + 2) / 3) + 1 )
            initial_index = int(svec + 8)
            threeatomlist = atom_coord_list[atom_index - 2: atom_index + 1]
            s_vec = angle_S_vector(threeatomlist)
            if svec - num_s_vectors[-3] == 0:
                array[initial_index - 9:] = s_vec[:]
            else:
                array[initial_index - 9: initial_index] = s_vec[:]
            s_vectors.append(array)
        if (svec - 2) % 3 == 0:
            #bonds
            array = np.zeros((len(atom_coord_list) * 3))
            atom_index = int(((svec + 1) / 3) + 2)
            initial_index = int(svec + 7)
            twoatomlist = atom_coord_list[atom_index - 2: atom_index]
            s_vec = length_S_vector(twoatomlist)
            if svec - num_s_vectors[-2] == 0:
                array[initial_index - 6:] = s_vec[:]
            else:
                array[initial_index - 6: initial_index] = s_vec[:]
            s_vectors.append(array)
        if svec % 3 == 0:
            #dihedrals
            array = np.zeros((len(atom_coord_list) * 3))
            atom_index = int((svec / 3) + 1)
            initial_index = int(svec + 6)
            fouratomlist = atom_coord_list[atom_index - 3: atom_index + 1]
            s_vec = dihedral_S_vector(fouratomlist)
            if svec - num_s_vectors[-1] == 0:
                array[initial_index - 12:] = s_vec[:]
            else:
                array[initial_index - 12: initial_index] = s_vec[:]
            s_vectors.append(array)
    bmat = np.stack(s_vectors[:])
    print("B Matrix successfully calculated")
    return bmat
