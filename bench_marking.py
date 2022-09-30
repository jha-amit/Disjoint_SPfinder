# -*- coding: utf-8 -*-
"""
Created on Mon Dec 20 18:15:03 2021

@author: Amit
"""
import math
import random
import numpy as np
import sys
import numpy as np
import random
import math
import csv
from cython_Dijkstra import node_val
from Cost_modified import interpolation,Cost_diamondgraph,Patchpoints,hash1,complete_graph,Patch_terminal,Cost_radial,interpolation_radial
import pickle
from bellman_ford import bellman_ford
from IDW import IDW
from matplotlib import pyplot as plt
import json
import time
from mpl_toolkits.mplot3d import Axes3D

#start[0]=57.26277; start[1]=-1115425; end[0]=58.2252; end[1]=-111.4611; distance=107km; k=31 (D/2tan(alpha) +1), K=60, r=53.5km

def shortest_path1(start_lat,start_long,end_lat,end_long,K,nTimes,alpha,L,r,input_points,center):
    
    

    conversion=1/111 #conversion from change in latitude to km.
                    #conversion for km to latitude = (111.32*math.cos(lat[k]*math.pi/180)) 
    end=[end_lat,end_long]   
    start=[start_lat,start_long]

    if ((K%2!=0) or (K < 0)):
        print('enter even K')
    else:
        # This parameter is to control the actual edge length. it divides 1 unit in l parts so we
                #can shrink or expand the grid for our uses. here i am using it for optimal size of grid to acheiv
                #  fairly accurate cost interpolations. 
        N=int(K/2)+1
        Patchx=np.zeros((N,N))
        Patchy=np.zeros((N,N))
        Patchz=np.zeros((N,N))
    
        theta=math.atan((((end_lat-start_lat)*111)/((end_long-start_long)*111.2))) # need correction for eliptical co-ordinates
        # We need to provide user data as (Lat,Long,Cost) This user data will be in the form of scattered data
        # point. But we need the costs along the edges of diamond frame. I will interpolate 
        # 'Inverse distance weight' method the user input for getting costs at control points which is nodes
        #  of diamond graph.
        #print(theta)
        start_time=time.time()
        Patchx,Patchy= Patchpoints(N,theta,alpha,1)
        end_time=time.time()
        print(end_time-start_time)
        k=0
        start_time=time.time()
        Lat=[1]*N**2
        Long=[1]*N**2
        for i in range(0,N):
            for j in range(0,N):
                Lat[k]=Patchy[i,j]*conversion+start[0]       
                Long[k]=Patchx[i,j]/(111.32*math.cos(Lat[k]*3.14/180))+start[1]
                k=k+1
        end_time=time.time()
        print('Patch to lat lon', end_time-start_time)
    
# generating lat long for benchmarking
        
        radius_meters=r
        radius_degrees = radius_meters*2*math.pi / 111300
        
        random.seed(10)
        w = [radius_meters*(random.random()) for i in range(input_points)]
        t_forlong =  [math.cos(2*math.pi*(random.random())) for i in range(input_points)]
        t_forlat=[math.sin(2*math.pi*(random.random())) for i in range(input_points)]
        
        temp_long_off=[]
        lat_off=[]


        for number1, number2,number3 in zip(w, t_forlong,t_forlat):
            lat_off.append(center[0]+(number1*number3)/111100)
            temp_long_off.append(center[1]+ (number1 * number2)/(111320*math.cos(lat_off[-1]*3.14/180)))
            

       # print(temp_long_off)

        Longitude_input=temp_long_off
        Latitude_input=lat_off
        random.seed(11)
        Cost_input=[random.random() for i in range(input_points)]
        with open('C:/Users/Amit/Desktop/Working_test/application_shortest_path/shortest_path/media/Cost_matrix.csv', "w") as csv_file:
            writer = csv.writer(csv_file, delimiter=',',lineterminator='\n')
            writer.writerow(['Longitude','Latitude','Cost'])
            for i in range(0, len(Latitude_input)):
                writer.writerow([Longitude_input[i],Latitude_input[i],Cost_input[i]])
        #print(Cost_input)
        
        start_time=time.time()
        Latitude_input=np.array(Latitude_input)
        Longitude_input=np.array(Longitude_input)
        Cost_input=np.array(Cost_input)
        Lat=np.array(Lat)
        Long=np.array(Long)
        
        start_time=time.time()
        

        Cost=IDW(Latitude_input,Longitude_input,Cost_input,Lat,Long,.1,.1)
        end_time=time.time()
        print('IDW', end_time-start_time)
        #print(Cost,N)
        start_time=time.time()
        k=0
        for i in range(0,N):
            for j in range(0,N):
                Patchz[i,j]=Cost[k]
                k=k+1        
        r=(nTimes-1)*(N-1)+1

        #print(Patchz)
        #print(Cost_input)

        zvals=np.zeros((r,r))
        zvals=interpolation(nTimes,N,Patchz)
        end_time=time.time()        
        print('Interpolation', end_time-start_time)
        start_time=time.time()
        L_Diag=2*L*math.sin(alpha)
        Cost_horizontal,Cost_vertical=Cost_diamondgraph(N,zvals,nTimes,L,L_Diag)
        end_time=time.time()
        print('Cost matrix', end_time-start_time)      

        n=int(1/48*((K+2)*(K+4)*(2*K+6)))
        Hash=np.zeros((n,3))
        #Node co-ordinates 
        Hash=hash1(n,K)

        N_K=np.array([1.0]*(K+1)) # number of nodes in each layer Note that there are K layers
                                                            #between K+1 planes.
        i=1
        #Number of nodes in a layer
        start_time=time.time()
        while i<=K+1:
            if i<=K/2:
                N_K[i-1]=int(i*(i+1)/2)
            else:
                N_K[i-1]=int((K-i+2)*(K-i+3)/2)
            i=i+1
        d=2
        Node_val=node_val(Hash,Cost_horizontal,Cost_vertical,N_K,K,d)

        #print(end_time-start_time)
        sPath=np.zeros((K,3))
        sPath[K-1,:]=Node_val[n-1,:]
        for m in range(K-2,-1,-1):
            sPath[m,:]=Node_val[int(sPath[m+1,1]),:]
               
        #print(sPath)
        end_time=time.time()
        print('shortest path', end_time-start_time)
    fig = plt.figure(num=1, clear=True)
    ax = fig.gca(projection='3d')
    
    return sPath, Hash, Patchx, Patchy, conversion,start,K, theta, nTimes,d
#start[0]=57.26277; start[1]=-111.5425; end[0]=58.2252; end[1]=-111.4611; distance=107km; k=31 (D/2tan(alpha) +1), K=60, r=53.5km
#center=[58.40333333333333,-111.75333333333333]
# center_200=[57.857813,-113.31986055555]
center_50=[(57.024166+57.9055)/2,(-111.53611-111.185)/2]
# center_100=[57.44652744445,-112.40111055554999]
# center_400=[58.6169441111,-115.3158327778]
start_lat=57.024166
start_long=-111.53611
end_lat=57.9055
end_long=-111.185
K=10
nTimes=13
alpha=0.785
L=1
r=50000
input_points=10000
center=center_50

#end=[57.4416666667,-112.3869444444;57.8688888889,-113.2661111111;58.6941666667,-115.1036111111;60.2097222222,-119.0955555556]
sPath,Hash,Patchx,Patchy,conversion,start,K,theta,nTimes,d= shortest_path1(start_lat,start_long,end_lat,end_long,K,nTimes,alpha,L,r,input_points,center)
#shortest_path1(start_lat,start_long,end_lat,end_long,K,nTimes,alpha,L,r,input_points,center)::




start_node_path1=1

highvalue=10000
target_node_path1=48
N_start_path1=2500
B_path1=np.random.rand(2500,2500)
start_time=time.time()

sPath_terminal_start_path1 = bellman_ford(B_path1,N_start_path1,start_node_path1,target_node_path1,highvalue)
# sPath_terminal_start_path2 = bellman_ford(N_start_path1,start_node_path2,target_node_path1,B_path2,highvalue)
end_time=time.time()
print("SSP execution time: ", end_time-start_time)



   