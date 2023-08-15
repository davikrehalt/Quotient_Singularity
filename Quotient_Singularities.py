import numpy as np
import itertools as itl
from scipy.linalg import null_space,block_diag

def dot_prod(n,a,b):
    return sum([a[i]*b[i] for i in range(n)])
    
def grid(dimension, group_size, weight_list):
    #returns the dual lattice points within the cube between the origin
    #and (group_size,group_size,...,group_size) except for the origin
    outputlist=[]
    for i in range((group_size+1)**dimension)[1:]:
        coordinates = tuple([(i//(group_size+1)**d)%(group_size+1) for d in range(dimension)])
        if dot_prod(dimension,coordinates, weight_list)%group_size==0:
            outputlist.append(coordinates)
    return outputlist

def sum_tuples(a,b):
    return tuple([a[i]+b[i] for i in range(len(a))])

def trim_generators(pre_generators):
    previous_step=sorted(pre_generators)
    current_step=set(pre_generators)
    minkowski_double=[]
    counter=0
    while(counter<len(previous_step)):
        for j in previous_step[counter:]:
            current_step.discard(sum_tuples(previous_step[counter],j))
        previous_step=sorted(current_step)
        counter+=1

    return previous_step
    
def getT1(R,generators,test=False):
    N=len(generators)
    dimension=len(generators[0])
    E_coordinate_diag = np.zeros((0,0))
    E_coordinates = np.zeros((dimension,0))
    E1=dict()
    index_lookup = dict()
    if test:
        print('E1_test')
        print(N)
    for d in range(dimension):
        E1[d]=sorted(filter(lambda x:generators[x][d]<R[d],list(range(N))))
        if test:
            print(str(E1[d]) +' for '+str(d))
        if E1[d]:
            counter=0
            for point in E1[d]:
                #looks up what is the index of the point (referred to by its index in the list of generators) in the matrix E_coordinate_matrix
                index_lookup[(d,point)]=E_coordinate_diag.shape[1]+counter
                counter+=1
            current_inclusion_matrix=np.array([list(generators[i]) for i in E1[d]]).transpose()
            E_coordinate_diag=block_diag(E_coordinate_diag,current_inclusion_matrix)
            E_coordinates=np.hstack((E_coordinates,current_inclusion_matrix))
    assert E_coordinate_diag.shape[1] == E_coordinates.shape[1] #number of Ei's (summed)
    E2=dict()
    Second_differential = np.zeros((E_coordinate_diag.shape[1],0))
    if test:
        print('E2_test')
    for pair in itl.combinations(range(dimension),2):
        E2[pair]=sorted(set(E1[pair[0]]).intersection(set(E1[pair[1]])))
        if test:
            print(str(E2[pair]) +' for '+str(pair))
        if E2[pair]:
            if test:
                print(np.array([[int(j == index_lookup[(pair[0],i)]) - int(j == index_lookup[(pair[1],i)])
                                 for j in range(E_coordinate_diag.shape[1])] for i in E2[pair]]).transpose())
            Second_differential=np.hstack((Second_differential,np.array([[int(j == index_lookup[(pair[0],i)]) - int(j == index_lookup[(pair[1],i)])
                                                                            for j in range(E_coordinate_diag.shape[1])] for i in E2[pair]]).transpose()))

    if E_coordinates.size == 0:
        return 0
    KERNEL_DIM = E_coordinates.shape[1]-np.linalg.matrix_rank(E_coordinates.astype(float)) #shape 1 is number of Ei's summed
    Null_space_of_E_coordinates = null_space(E_coordinate_diag)

    Combined_matrix=np.hstack((Null_space_of_E_coordinates,Second_differential))

    if Combined_matrix.size == 0:
        IMAGE_DIM = 0
    else:
        IMAGE_DIM = np.linalg.matrix_rank(Combined_matrix)
    if test:
        print(KERNEL_DIM)
        print(IMAGE_DIM)
    return KERNEL_DIM-IMAGE_DIM

def compute(dimension, group_size, weight_list):
    Grid=grid(dimension, group_size, weight_list)
    generators=trim_generators(Grid)
    print('gen: '+str(generators))
    output=dict()
    zerotuple=tuple([0 for _ in range(dimension)])
    l=len(Grid)
    for (ind,R) in enumerate(Grid):
        if ind%100==0:
            print('computed: '+str(ind) + ' out of' + str(l))
        for subsetsize in range(dimension):
            for subset in itl.combinations(range(dimension),subsetsize):
                test_R = tuple([(-1 if i in subset else R[i]) for i in range(dimension)])
                T1_of_R=getT1(test_R,generators)
                if T1_of_R!=0:
                    output[test_R]=T1_of_R
    return output

if __name__ == "__main__":
    print(compute(3,121,(1,22,32)))
    

