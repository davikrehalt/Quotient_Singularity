import numpy as np
import itertools as itl
from scipy.linalg import null_space,block_diag
from sympy.ntheory import factorint
from sympy.core.numbers import igcd
from copy import copy

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
#Todo check if doing twice is worth it!
    
def getT1(R,generators,test=False):
    N=len(generators)
    dimension=len(generators[0])
    E_coordinate_diag = np.zeros((0,0))
    E_coordinates = np.zeros((dimension,0))
    E1=dict()
    index_lookup = dict()
    for d in range(dimension):
        E1[d]=sorted(filter(lambda x:generators[x][d]<R[d],list(range(N))))
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
    for pair in itl.combinations(range(dimension),2):
        E2[pair]=sorted(set(E1[pair[0]]).intersection(set(E1[pair[1]])))
        if E2[pair]:
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
    return KERNEL_DIM-IMAGE_DIM

def compute2d(group_size, weight_list):
    if group_size==1:
        return dict()
    gcd_list=[igcd(weight_list[i],group_size) for i in range(2)]
    Grid=grid(2, group_size, weight_list)
    generators=trim_generators(Grid)
    output=dict()
    l=len(Grid)
    for (ind,R) in enumerate(Grid):
        T1_of_R=getT1(R,generators)
        if T1_of_R!=0:
            output[R]=T1_of_R
    return output

def check_interior(dimension, group_size, weight_list):
    if dimension != 3:
        return NotImplementedError
    #pre-check
    if igcd(weight_list[0],weight_list[1],group_size)!=1:
        raise ValueError
    if igcd(weight_list[0],weight_list[2],group_size)!=1:
        raise ValueError
    if igcd(weight_list[1],weight_list[2],group_size)!=1:
        raise ValueError
    #check if more than two prime factors
    if group_size<6:
        return dict()
    factors=factorint(group_size)
    if len(factorint(group_size))<2:
        return dict()
    gcd_list=[igcd(weight_list[i],group_size) for i in range(3)]
    if sum([int(gcd_list[i]>1) for i in range(3)])<2:
        return dict()

    #generate grid
    Grid=[]
    for i in range(group_size//gcd_list[0]+1):
        for j in range(group_size//gcd_list[1]+1):
            for k in range(group_size//gcd_list[2]+1):
                if (weight_list[0]*i+weight_list[1]*j+weight_list[2]*k)%group_size==0:
                    if (i,j,k)!=(0,0,0):
                        Grid.append((i,j,k))
    generators=trim_generators(Grid)

    R_dict=dict()
    l=len(Grid)
    for (ind,R) in enumerate(Grid):
        if R[0]>0 and R[1]>0 and R[2]>0:
            #all positive
            tocheck = set()
            if R[0]>1 and gcd_list[1]>=R[0] and gcd_list[2]>=R[0]:
                tocheck.add(0)
            if R[1]>1 and gcd_list[0]>=R[1] and gcd_list[2]>=R[1]:
                tocheck.add(1)
            if R[2]>1 and gcd_list[0]>=R[2] and gcd_list[1]>=R[2]:
                tocheck.add(2)
            if len(tocheck)>0:
                if len(tocheck)>0:
                    T1_of_R=getT1(R,generators)
                    if T1_of_R>0:
                        R_dict[R]=T1_of_R
    return R_dict

def compute3d(group_size,weight_list):
    assert group_size>1
    assert len(weight_list)==3
    assert weight_list[0]>=0
    assert weight_list[1]>=0
    assert weight_list[2]>=0
    #fast threefold calculation
    Interior = check_interior(3,group_size,weight_list) #interior R's

    gcd_list=[igcd(weight_list[i],group_size) for i in range(3)]
    Boundaries=[]
    #Now we compute what the boundaries are
    weights = (weight_list[1]%gcd_list[0],weight_list[2]%gcd_list[0])
    Boundaries.extend(map(lambda x: ((-1,x[0][0],x[0][1]),x[1]),list(compute2d(gcd_list[0],weights).items())))
    weights = (weight_list[0]%gcd_list[1],weight_list[2]%gcd_list[1])
    Boundaries.extend(map(lambda x: ((x[0][0],-1,x[0][1]),x[1]),list(compute2d(gcd_list[1],weights).items())))
    weights = (weight_list[0]%gcd_list[2],weight_list[1]%gcd_list[2])
    Boundaries.extend(map(lambda x: ((x[0][0],x[0][1],-1),x[1]),list(compute2d(gcd_list[2],weights).items())))
    return dict(list(Interior.items())+Boundaries)

if __name__ == "__main__":
    print(compute3d(99,(1,22,90)))
    

