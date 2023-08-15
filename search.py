from Quotient_Singularities import *
from sympy.ntheory import factorint
from sympy.core.numbers import igcd
from copy import copy

def has_positives(dimension, group_size, weight_list):
    if dimension != 3:
        return NotImplementedError
    #pre-check
    if igcd(weight_list[0],weight_list[1],group_size)!=1:
        return False
    if igcd(weight_list[0],weight_list[2],group_size)!=1:
        return False
    if igcd(weight_list[1],weight_list[2],group_size)!=1:
        return False
    #check if more than two prime factors
    if group_size<6:
        return False
    factors=factorint(group_size)
    if len(factorint(group_size))<2:
        return False
    gcd_list=[igcd(weight_list[i],group_size) for i in range(3)]
    if sum([int(gcd_list[i]>1) for i in range(3)])<2:
        return False

    #generate grid
    Grid=grid(dimension, group_size, weight_list)
    generators=trim_generators(Grid)
    R_list=dict()
    zerotuple=tuple([0 for _ in range(dimension)])
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
                for g in generators:
                    if len(tocheck)==0:
                        break
                    new_tocheck=copy(tocheck)
                    for i in tocheck:
                        if g[i]<R[i] and g[i]>0:
                            if sum([int(g[i]<R[i]) for i in range(3)])>1:
                                new_tocheck.remove(i)
                    tocheck=new_tocheck
                if len(tocheck)>0:
                    return R
    return False

if __name__ == "__main__":
    for group_size in range(100)[6:]:
        print('now calculating n='+str(group_size))
        factors=list(factorint(group_size))
        if len(factorint(group_size))>=2:
            for A in range(group_size)[2:]:
                for B in range(group_size)[A:]:
                    for C in range(group_size)[B:]:
                        out=has_positives(3,group_size,(A,B,C))
                        if out:
                            print('1/'+str(group_size)+str((A,B,C)))
                            print('has R')
                            print(out)
                    
            
    
