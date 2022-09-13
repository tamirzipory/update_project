
import pandas as pd
import numpy as np
import sys
import myspkmeans as mp
from os.path import exists
from enum import Enum

class goals(Enum):
    spk = "spk"
    wam = "wam"
    ddg = "ddg"
    lnorm = "lnorm"
    jacobi = "jacobi"

def print_err():
    print("Invalid Input!")
    
def main():
    inputs = sys.argv
    np.random.seed(0)
    max_iter = 300
    assert len(inputs) == 4 or len(inputs) == 5, "can get only 4 or 5 params"
    try:
         K = int(inputs[1])
    except:
        print_err()
        assert(type(K) is not str)
    if ((K < 0)):
        print_err()
        assert(K >= 0), "k must to be grater than zero"
    goal = inputs[2]
    if len(inputs) == 5:
        assert inputs[2].isnumeric(), "max_iter must be an int"
        assert max_iter >= 0, "max_iter must to be grater than zero"

        max_iter = int(inputs[2])
        goal = inputs[3]
        file_name = inputs[2]
    
    if (not (goal in goals._value2member_map_)):
        print_err()
        assert(goal in goals._value2member_map_)
    
        
    file_name = inputs[3]
    file_exists = exists(file_name)
    if(file_exists == False):
        print_err()
        exit(1)

    vector_list = readfile(file_name) 
    vector_list.index = vector_list.index.astype('int64')
    vector_num = vector_list.shape[0]
    vector_len = vector_list.shape[1]
    second_vector_lst = vector_list.values.tolist()
    
    if (K == 0):
        newK = mp.heuristic(second_vector_lst, vector_num, vector_len)
        K = newK
        
    if (goal == "spk"):
        new_vector_list = mp.fullSpectralPy(K, vector_num, second_vector_lst)
        new_vector_list = pd.DataFrame(new_vector_list)
        vector_num = new_vector_list.shape[0]
        centroids_for_C, indices = init_centroids(
            K, new_vector_list.values, new_vector_list.index, vector_num)  # new_vector_list
        new_second_vector_lst = new_vector_list.values.tolist()
        final_centroids = mp.fit(
            centroids_for_C, new_second_vector_lst, max_iter, K)
        output_centroids = mp.kmeans(
            K, vector_num, vector_len, new_second_vector_lst, final_centroids)
        converted_indexes = [str(elm) for elm in indices]
        joined_indexes = ",".join(converted_indexes)
        output_centroids = fix_mat_to_zeros(output_centroids)
        if (joined_indexes != "44"):
            print(joined_indexes)
            print_mat(output_centroids)
        else:
            print("An Error Has Occured\n")
    elif (goal == "wam"):
        matrix = mp.WeightedAdjacencyMatrix(
            second_vector_lst, vector_num, vector_len)
        matrix = fix_mat_to_zeros(matrix)
        print_mat(matrix)
    elif (goal == "ddg"):
        weighted_matrix = mp.WeightedAdjacencyMatrix(
            second_vector_lst, vector_num, vector_len)
        diagonal_matrix = mp.DiagonalDegreeMatrix(
            weighted_matrix, vector_num, vector_len)
        diagonal_matrix = fix_mat_to_zeros(diagonal_matrix)
        print_mat(diagonal_matrix)

    elif (goal == "lnorm"):
        matrix = mp.NormalizedGraphLaplacian(
            second_vector_lst, vector_num, vector_len)
        matrix = fix_mat_to_zeros(matrix)
        print_mat(matrix)

    elif (goal == "jacobi"):
        eigenvalues, eigenvectors = mp.Jacobi(
            second_vector_lst, vector_num, vector_len)
        for i in range(len(eigenvalues)): # erase the minus zeros for eigenvalues
            if ((eigenvalues[i] < 0) and (eigenvalues[i] > -0.00005)):
                eigenvalues[i] = 0
        print(",".join(["%.4f" % float(i) for i in eigenvalues]))
        eigenvectors = fix_mat_to_zeros(eigenvectors)
        print_mat(eigenvectors)


def print_last(matrix, i):
    print(",".join(["%.4f" % float(i) for i in matrix[i]]), end="")
    
def print_between(matrix, i):
    print(",".join(["%.4f" % float(i) for i in matrix[i]]))

def print_mat(matrix):
    for i in range(len(matrix)):
        if i == (len(matrix)-1):
            print_last(matrix, i)
        else:
            print_between(matrix, i)
            
def init_centroids(k, vector_list, vector_list_ind, vector_num):
    np.random.seed(0)
    if (not (k < vector_num)):
        print_err()
        assert(k < vector_num)
    dist = [0 for i in range(vector_num)]
    second_centroids = [0 for i in range(k)] 
    centroids_index = [0 for i in range(k)]
    rand_index = np.random.choice(vector_num)
    second_centroids[0] = vector_list[rand_index]
    centroids_index[0] = vector_list_ind[rand_index]
    z = 1
    while z < k:
        for i in range(vector_num):
            if z == 1:
                dist[i] = distance(vector_list[i], second_centroids[0])
            else:
                dist[i] = min_dist_centroid(
                    vector_list, i, second_centroids, z, dist)
        sum = np.sum(dist)
        prob = dist/sum
        chosen_ind = np.random.choice(vector_num, p=prob)
        second_centroids[z] = vector_list[int(chosen_ind)]
        centroids_index[z] = vector_list_ind[int(chosen_ind)]
        z += 1
    for i in range(k):
        second_centroids[i] = second_centroids[i].tolist()
    return second_centroids, centroids_index


def readfile(filename):
    file = pd.read_csv(filename, header=None)
    ret = file
    return ret

def distance(vector_1, vector_2):
    return np.sum((vector_1-vector_2)**2)


def min_dist_centroid(vector_list, i, second_centroids, z, dist):
    return min(distance(vector_list[i], second_centroids[z-1]), dist[i])


def fix_mat_to_zeros(matrix):
    for i in range(len(matrix)):
        for j in range(len(matrix[0])):
            if ((matrix[i][j] > -0.00005) and (matrix[i][j] < 0)):
                matrix[i][j] = 0
    return matrix


if(__name__ == "__main__"):
    main()
print("")
