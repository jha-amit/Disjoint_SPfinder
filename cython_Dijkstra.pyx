# -*- coding: utf-8 -*-
"""
Created on Thu Feb 11 16:24:53 2021

@author: amit
"""


# -*- coding: utf-8 -*-
"""
Created on Mon Sep 14 16:46:18 2020

@author: amit
"""
import sys
import numpy as np
import random
import math
cimport cython
import time
 
ctypedef fused my_type:
    int
    double
    #long long

#to save us some more time as we are sure about indices not going out of
#bounds and there are no negative index to deal with so we can supress these
#two processes 
@cython.boundscheck(False) # Deactivate bounds checking
@cython.wraparound(False) # Deactivate negative indexing.

#Here we are type declaring the Hash and other variables using fused type function
def node_val(my_type[:, ::1] Hash, my_type[:, ::1] Cost_horizontal, my_type[:, ::1] Cost_vertical, my_type[:, ::1] Cost_diag, my_type[::1] N, my_type[::1] truncated_layer_start_cost,my_type[::1] truncated_layer_end_cost,my_type K, my_type d, truncated_layer):
    
    n = Hash.shape[0]
    k = Hash.shape[1]
    N_nodes=int(K/2) +1
    #assert tuple(.shape) == tuple(array_2.shape)
#Here accprding to 
    
    #elif my_type is cython.longlong:
        #dtype = np.longlong
    
    if my_type is int:
        dtype = np.intc
    elif my_type is double:
        dtype = np.double
        
    Node_val1 = np.random.rand(n, k)*np.inf
    #Node_val2 = np.random.rand(n, k)*np.inf
    
    Node_val1[0,0]=0.0
    Node_val1[0,1]=0.0
    Node_val1[0,2]=0.0
    
    
    cdef double[:,::1] node_val1 = Node_val1
    #cdef double[:,::1] node_val2 = Node_val2
   
    cdef double  ip, jp, i1, i2, j, l1, l2
    cdef double Cijk = 0, Cijk1 = 0
    cdef int v1, v2,j1,j2,j3,source_neighbour,sink_neighbour,i,layer

    source_neighbour = truncated_layer
    sink_neighbour = int(K-source_neighbour)

    Cijk=0
    i=0
    
    for v2 in range(int(sum(N[:source_neighbour])),int(sum(N[:source_neighbour+1]))):
        i1 = Hash[v2,1]
        i2 = Hash[v2,2]
       
        if abs(i1-i2)>=d:
            node_val1[v2,0] = truncated_layer_start_cost[i]+truncated_layer_start_cost[int(i+(i2-i1)/2)] # i have to read the indices reverse to sync with the Hash (3D graph nodes)
            node_val1[v2,1] = 0
            node_val1[v2,2] = v2
            #print(node_val1[v2,0],v2)
            i = i+1
          

    layer = source_neighbour

    while layer <=sink_neighbour:
        
        res = all(ele == math.inf for ele in node_val1[int(sum(N[:layer])):int(sum(N[:layer+1])),0])
        
        if res == True:
             print('shift the target, its in infinite cost zone or give a bigger area for IDW point interpolation.')
             break
         
        if layer < sink_neighbour:   
            for v1 in range(int(sum(N[:layer])),int(sum(N[:layer+1]))):        
                # give range of nodes in next layer.
                j=(Hash[v1,0])            
                if layer < sink_neighbour-1:
                    l=int(sum(N[:layer+3]))
                else:
                    l=int(sum(N[:layer+2]))
                    
                
                for v2 in range(int(sum(N[:layer+1])),l):
                    #print(l1,l2)
                    j=Hash[v2,0]
                    
                    #is there an edge from v1 to v2, namely
                    #determine if an edge exist to j-level vertex v2
                     #from the previous(j-1)-level vertex v1
                      #at least the vertices are at the right j-levels
                        #get the j coordinate of vertex v2
                        #potentially there is an edge, but need to inevstigate further
                        #namely, check if it is
                        #too early/late to differentiate the paths OR
                        #the vertices are far apart (i-distance is at least 2)
                    if ((abs(Hash[v2,1]-Hash[v2,2]) >= d) and\
                        (abs(Hash[v2,1]-Hash[v1,1]) <= 1) and\
                        (abs(Hash[v2,2]-Hash[v1,2]) <=1)):
                        #in this case, add the edge
                        #determine its (cumulative) cost first                        
                        i1 = (Hash[v2,1])
                        i2 = (Hash[v2,2])

                        ip = 0.5*(j-i1)
                        jp = 0.5*(j+i1)

                        if abs(Hash[v2,0]-Hash[v1,0])<=1:

                            if Hash[v1,1]<=Hash[v2,1]: 
                                Cijk = Cost_horizontal[int(ip),int(jp)]
                            else:
                                Cijk = Cost_vertical[int(ip),int(jp)]

                        elif abs(Hash[v2,0]-Hash[v1,0])<=2:

                            Cijk = Cost_diag[int(ip),int(jp)]

                            if (Cijk >= Cost_horizontal[int(ip),int(jp)] + Cost_vertical[int(ip),int(jp) - 1])\
                                        and (Cijk >= Cost_vertical[int(ip),int(jp)] + Cost_horizontal[int(ip)-1,int(jp)]):

                                Cijk = np.inf
                                #max((Cost_horizontal[int(ip),int(jp)] + Cost_vertical[int(ip),int(jp) - 1]),\
                                        #(Cost_vertical[int(ip),int(jp)] + Cost_horizontal[int(ip)-1,int(jp)]))


                        ip = 0.5*(j-i2)
                        jp = 0.5*(j+i2)

                        if abs(Hash[v2,0]-Hash[v1,0])<=1:                                             

                            if Hash[v1,2]<=Hash[v2,2]: 
                                Cijk = Cost_horizontal[int(ip),int(jp)]+Cijk # changing the inquality here
                            else:
                                Cijk = Cost_vertical[int(ip),int(jp)]+Cijk

                        elif abs(Hash[v2,0]-Hash[v1,0])<=2:

                            Cijk1 = Cost_diag[int(ip),int(jp)]

                            if (Cijk1 >= Cost_horizontal[int(ip),int(jp)] + Cost_vertical[int(ip),int(jp) - 1])\
                                        and (Cijk1 >= Cost_vertical[int(ip),int(jp)] + Cost_horizontal[int(ip)-1,int(jp)]):

                                Cijk1 = np.inf

                            Cijk = Cijk1+Cijk


                        if (node_val1[v2,0]>(Cijk+node_val1[v1,0])):
                    
                            node_val1[v2,0]=Cijk+node_val1[v1,0]                            
                            node_val1[v2,1]=v1
                            node_val1[v2,2]=v2

        else:
            i=0
            for v1 in range(int(sum(N[:layer])),int(sum(N[:layer+1]))):
                i1=Hash[v1,1]
                i2=Hash[v1,2]
                if abs(i1-i2)>=d:

                    Cost_truncated_end1 = truncated_layer_end_cost[i]+ truncated_layer_end_cost[int(i+(i1-i2)/2)]
                    
                    if (node_val1[n-1,0]>(Cost_truncated_end1 + node_val1[v1,0])):           
                        node_val1[int(n-1),0]=Cost_truncated_end1+node_val1[v1,0]
                        node_val1[int(n-1),1]=v1
                        node_val1[int(n-1),2]=float(n-1)
                        i=i+1                       
                                
        layer = layer+1
    
    return Node_val1
                


#if ((v1//N + v1%N) >= (v1//N + d/2) and (v1//N <= v1%N)) or ((v1//N + v1%N) \
#                                                                    <= (v1//N + d/2) and (v1//N >= v1%N)): 
    