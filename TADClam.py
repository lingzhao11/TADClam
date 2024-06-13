"""TADClam module"""

import numpy as np
import overlap_multitree as MT
from multiprocessing import Pool, Manager
from scipy.special import digamma
import argparse

def sigm(x):
    """
    digamma function of vlue x
    """
    if x > 0:
        return  np.divide(np.exp(-1. * x), 1. - np.exp(-1. * x))
    else:
        return 0

def gradient(F, A, i):
    """
    calculating the gradient of table F
    """
    N, C = F.shape
    #neighbor sum/3
    sum_neigh = np.zeros((C,))
    for nb in range(i- int(N/ 3), i+ int(N/ 3), 1):
        if nb != i and nb >= 0 and nb < N:
            dotproduct = F[nb].dot(F[i]) + 1
            fi = digamma(dotproduct)
            sum_neigh += F[nb] * (np.log(A[i, nb]) - fi)  

    return sum_neigh


def train_multipocess(A, C, iterations, my_list):
    """
    process for calculating table F
    """
    N = A.shape[0]
    np.random.seed(0)
    F = np.random.uniform(0, 1, (N, C))
    for n in range(iterations):
        for person in range(N):
            grad = gradient(F, A, person)
            F[person] += 0.006 * grad
            F[person] = np.maximum(0.001, F[person])     
    # if < mean(f), =0    
    F_mean= np.mean(F)
    for i in range(len(F)):
        for j in range(C):
            if F[i][j] < F_mean:
                F[i][j] = 0
    # select the first and second score
    for i in range(len(F)):
        tmp_list = sorted(F[i])
        for j in range(C):
            if F[i][j] < tmp_list[-2]:
                F[i][j] = 0 

    temp_boundary = get_boundary(F, C)

    for new_TAD in temp_boundary:
        if new_TAD not in my_list:
            my_list.append(new_TAD)

def first_and_last_index(li, lower_limit = 0):
    """
    selecting candidate TADs
    """
    result = []
    foundstart = False
    foundend = False
    startindex = 0
    endindex = 0
    for i in range(0, len(li)):
        if li[i] > lower_limit :
            if not foundstart:
                foundstart = True
                startindex = i
        else:
            if foundstart:
                foundend = True
                endindex = i - 1
        if foundend:
            result.append([startindex, endindex])
            foundstart = False
            foundend = False
            startindex = 0
            endindex = 0
    if foundstart:
        result.append([startindex, len(li)-1])
    return result


def get_boundary(F, k): 
    #return TADs >= 3 bins
    boundaries= []
    for cluster in range(k):
        temp_arr = F[ :, cluster]
        boundary_list = first_and_last_index(temp_arr)
        for temp_boundary in boundary_list:
            if temp_boundary[1]- temp_boundary[0] >= 3:
                boundaries.append([temp_boundary[0], temp_boundary[1]])         
    return boundaries



def candidate_find_multi_process(adj, iter, max_k):
    """
    multi-process for calculating table F
    """
    
    manager = Manager()
    my_list = manager.list()
    pool = Pool(processes = max_k - 2)
    for k in range(1, max_k + 1, 1):
        pool.apply_async(train_multipocess, (adj, k, iter, my_list,))
    pool.close()
    pool.join()
    return list(my_list)


def boundary_sorted(N, mid_boundaries):
    # for each bin, sorted the tad according to the len of tad, big to small
    sorted_boundaries = []               
    for i in range(N):
        parent_i = []
        len_parent_i = []
        for temp in mid_boundaries:
            if temp[0] <= i and  i <= temp[1]:
                parent_i.append(temp)
                len_parent_i.append(temp[1]-temp[0])
        sorted_indices = sorted(range(len(len_parent_i)), key=lambda x: len_parent_i[x], reverse=True)
        for index in sorted_indices:
            if not parent_i[index] in sorted_boundaries:
                sorted_boundaries.append(parent_i[index])
    return sorted_boundaries
      
def filtering_noise(adj, cadidate_boundaries, T_on = 2):
    #basic selecting, filting noise
    boundaries = cadidate_boundaries
    while True:
        new_cadidate_boundaries = boundaries.copy()
        temp_boundaries = []
        score = []
        for temp in boundaries:
            inter_mean = np.mean(adj[temp[0]: temp[1]+ 1,temp[0]: temp[1]+ 1])
            len_tad = temp[1]- temp[0] + 1
            #left part mean
            left_mean = 0
            over_left_mean= 0
            non_over_left_mean = 0
            if temp[0] > 0:
                if temp[0] - len_tad > 0:
                    left_start = temp[0]- len_tad
                else:
                    left_start = 0
                non_over_left_mean = np.mean(adj[temp[0]: temp[1]+ 1, left_start: temp[0]])
                #if overlap 
                for lap_tad in boundaries:
                    if lap_tad[0]< temp[0] and temp[0]< lap_tad[1] and lap_tad[1]< temp[1]:
                        over_left_mean = max(over_left_mean, np.mean(adj[lap_tad[1]+ 1: temp[1]+ 1, left_start: temp[0]]))

            if over_left_mean>0:
                left_mean= over_left_mean
            else:
                left_mean= non_over_left_mean
            #right part mean          
            right_mean = 0
            over_right_mean = 0
            non_over_right_mean = 0
            if temp[1] < len(adj)- 1:
                if temp[1]+ len_tad < len(adj):
                    right_end= temp[1]+ len_tad
                else:
                    right_end= len(adj) - 1
                non_over_right_mean = np.mean(adj[temp[0]: temp[1]+ 1, temp[1]+ 1: right_end+1])
                # if overlap 
                for lap_tad in boundaries:
                    if lap_tad[0] < temp[1] and temp[1]< lap_tad[1] and temp[0]< lap_tad[0]:
                        over_right_mean = max(over_right_mean,  np.mean(adj[temp[0]: lap_tad[0], temp[1]+ 1:right_end+1]))
            if over_right_mean>0:
                right_mean= over_right_mean
            else:
                right_mean= non_over_right_mean         
            
            left_right_mean = max(left_mean, right_mean)
            
            if left_right_mean > 0:
                if inter_mean/ left_right_mean > T_on:
                    if len(temp_boundaries) == 0:
                        temp_boundaries.append(temp)
                        score.append(inter_mean/ left_right_mean)
                    else:
                        sin = 0
                        #if two tad has the same simularity,choose the best
                        single = -1
                        for i in range(len(temp_boundaries)):
                            if abs(temp[0]-temp_boundaries[i][0]) <= 3 and  abs(temp[1]-temp_boundaries[i][1]) <= 3 and inter_mean/ left_right_mean > score[i]:
                                sin = 1
                                single = i
                                break
                            if abs(temp[0]- boundaries[i][0]) <= 2 and  abs(temp[1]- boundaries[i][1]) <= 2 and inter_mean/ left_right_mean <= score[i]:
                                sin = 2
                                break
                        # no simular tad
                        if sin == 0:
                            temp_boundaries.append(temp) 
                            score.append(inter_mean/ left_right_mean)
                        # has simular tad, but this one is best than the old
                        if sin == 2:
                            if temp in boundaries:
                                boundaries.remove(temp)
                        if sin == 1:
                            if temp_boundaries[single] in boundaries:
                                boundaries.remove(temp_boundaries[single])
                            temp_boundaries.remove(temp_boundaries[single])
                            score.remove(score[single])
                            temp_boundaries.append(temp) 
                            score.append(inter_mean/ left_right_mean)  
        if new_cadidate_boundaries == temp_boundaries:
            break
        else:
            boundaries = temp_boundaries


    return boundaries 


def score_filter(adj, mid_boundaries, T_div= 2.5):
    #score filtering
    sorted_boundaries = boundary_sorted(len(adj), mid_boundaries)
    Multitree = MT.MultiChildTree(0, len(adj) - 1)
    for temp in sorted_boundaries:
        Multitree.insert(temp[0], temp[1])
    node_list = Multitree.acquire_list()  
    node_list = list(set(node_list))
    final_boundaries = []
    while True:
        new_node_list = node_list.copy()
        for node in node_list:
            if node.parent != None:
               
                for node_parent in node.parent:
                    # inter _node
                    node_inter=[]
                    for i in range(node.val[0], node.val[1] + 1):
                        for j in range(node.val[0], node.val[1] +1):
                            point = [i,j]
                            node_inter.append(point)
                    node_inter_orign = node_inter
                    # child_inter_node
                    node_child_inter_sum=[]
                    for node_child in node.child:
                        node_child_inter=[]
                        for i in range(node_child.val[0], node_child.val[1] + 1):
                            for j in range(node_child.val[0], node_child.val[1] +1):
                                point = [i,j]
                                node_child_inter.append(point)
                        node_child_inter_sum = [val for val in node_child_inter if val not in node_child_inter_sum] + node_child_inter_sum
                    #overlap_node
                    overlap_node_list = find_overlap_node(node, node_list)
                    node_overlap_inter_sum=[]
                    for overlap_node in overlap_node_list:
                        node_overlap_inter=[]
                        for i in range(overlap_node.val[0], overlap_node.val[1] + 1):
                            for j in range(overlap_node.val[0], overlap_node.val[1] +1):
                                point = [i,j]
                                node_overlap_inter.append(point)
                        node_overlap_inter_sum = [val for val in node_overlap_inter if val not in node_overlap_inter_sum] + node_overlap_inter_sum
                    
                    node_inter = [val for val in node_inter if val not in node_child_inter_sum]
                    node_inter = [val for val in node_inter if val not in node_overlap_inter_sum]
                    #node average
                    node_value_sum = 0
                    for node_point in node_inter:
                        node_value_sum = node_value_sum + adj[node_point[0], node_point[1]]
                    if len(node_inter) > 0:
                        node_value_ave = node_value_sum / len(node_inter)
                    else:
                        node_value_ave = 0
                    # outer bin 
                    outer_bin = []
                    for node_point in node_inter: 
                        if node_point[0] not in outer_bin:
                            outer_bin.append(node_point[0])
                    
                    outer_point = []
                    for bin in outer_bin:
                        for j in range(node_parent.val[0], node_parent.val[1] +1): 
                            point = [bin, j]
                            outer_point.append(point)

                    outer_point = [val for val in outer_point if val not in node_overlap_inter_sum]
                    outer_point = [val for val in outer_point if val not in node_inter_orign]
                    #node outer average
                    node_outer_sum = 0
                    for node_point in outer_point:
                        node_outer_sum = node_outer_sum + adj[node_point[0], node_point[1]]
                    if len(outer_point) > 0:
                        outer_ave = node_outer_sum / len(outer_point)
                    else:
                        outer_ave = 0
                    if outer_ave > 0:
                        DIV =  node_value_ave / outer_ave
                    else:
                        DIV = 0
                    if DIV < T_div:
                        for i in range(0, len(node.child), 1):
                            node.child[i].parent.append(node_parent)
                            if node in node.child[i].parent:
                                node.child[i].parent.remove(node)
                            node_parent.child.append(node.child[i])
                        if node in node_parent.child:
                            node_parent.child.remove(node)
                        if node in node_list:
                           node_list.remove(node)

        if new_node_list == node_list:
            break
        else:
            mid_boundaries = []
            for node in node_list:
                if node.val[1]- node.val[0] >= 4:
                    mid_boundaries.append([node.val[0], node.val[1]])
            sorted_boundaries = boundary_sorted_v2(len(adj), mid_boundaries)
            Multitree = MT.MultiChildTree(0, len(adj) - 1)
            for temp in sorted_boundaries:
                Multitree.insert(temp[0], temp[1])
            node_list = Multitree.acquire_list()  
            node_list = list(set(node_list))

    for node in node_list:
        final_boundaries.append([node.val[0], node.val[1]])

    final_boundaries_v2 = []
    for temp in final_boundaries:
        x = [temp[0] + 1, temp[1] + 1]
        final_boundaries_v2.append(x)

    return  final_boundaries


def find_overlap_node(node, node_list):
    overlap_node_list = []  

    for temp in node_list:
        if (temp.val[0] < node.val[0] < temp.val[1] < node.val[1]) or ( node.val[0] < temp.val[0] < node.val[1] < temp.val[1]):
            overlap_node_list.append(temp)
    return overlap_node_list

def boundary_sorted_v2(N, mid_boundaries):
    # for each bin, sorted the tad according to the len of tad, big to small
    sorted_boundaries=[]               
    for i in range(N - 1, -1, -1):
        parent_i = []
        len_parent_i = []
        for temp in mid_boundaries:
            if temp[0] <= i and  i <= temp[1]:
                parent_i.append(temp)
                len_parent_i.append(temp[1]-temp[0])
        sorted_indices = sorted(range(len(len_parent_i)), key=lambda x: len_parent_i[x], reverse=True)
        for index in sorted_indices:
            if not parent_i[index] in sorted_boundaries:
                sorted_boundaries.append(parent_i[index])
    return sorted_boundaries




def run(matrix, iter = 1500, thr = 2.50):

    adj_original = np.loadtxt(matrix)
    max_k = int(len(adj_original)**0.5)
    eval_matrix = adj_original/ adj_original.max()
    min = 100000
    for i in range(len(adj_original)):
        for j in range(len(adj_original)):
            if adj_original[i][j]> 0:
                if adj_original[i][j]< min:
                    min = adj_original[i][j]
    adj = adj_original + min
    np.fill_diagonal(adj, 0)
    cadidate_boundaries = candidate_find_multi_process(adj, iter, max_k)
    mid_boundaries = filtering_noise(eval_matrix, cadidate_boundaries)
    tad = score_filter(eval_matrix, mid_boundaries, thr)
    return tad

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-in", dest="input", type=str, default="test.txt", help="destination of the data")
    parser.add_argument("-iter", dest="iter",type=int, default=1500, help="iteration of table F")
    parser.add_argument("-threshold", dest="Threshold", type=float, default=2.50, help="Threshold of filtering")
    parser.add_argument("-out", dest="out", type=str, default="result.txt", help="output of the results")
    args = parser.parse_args()
    TADs = run(args.input, args.iter, args.Threshold)
    np.savetxt(args.out,TADs)

if __name__ == "__main__":
    main()
