from django.shortcuts import redirect, render
import sys
import numpy as np
import random
import math
import csv
import os

from cython_Dijkstra import node_val
from Cost_modified import interpolation,Cost_diamondgraph,Patchpoints,hash1,Patch_terminal,interpolation_radial,Cost_radial
import pickle
from .forms import UploadFileForm
from django.conf import settings
from django.core.files.storage import FileSystemStorage
from django.core.files.storage import FileSystemStorage
from django.db import models
from numpy import genfromtxt
import json 
from django.core.files.storage import default_storage
from IDW import IDW
from bellman_ford import bellman_ford
from django.views.decorators.csrf import csrf_exempt
from django.http import HttpResponse, JsonResponse
from django.core import serializers


val=None
def distance(x1,y1,x2,y2):
        dist=math.sqrt((x1-x2)**2+(y1-y2)**2)
        return dist

def home(request):
   return render(request,'Home_test1.html',{"name":"user"}) 

def plot(request):
   
    if  request.POST:
        start_lat=float(request.POST['Start_point_latitude']) 
        start_long=float(request.POST['Start_point_longitude'])
        end_lat=float(request.POST['End_point_latitude'])  
        end_long=float(request.POST['End_point_longitude'])
        d=int(request.POST['distance'])
        alpha=float(request.POST['edge_slope'])
        gmap_license=request.POST["enter_your_gmap_API_key"]
        l_scale=float(request.POST["enter scaling factor for grid"])


    conversion=1/111 #conversion from change in latitude to km.
                    #conversion for km to latitude = (111.32*math.cos(lat[k]*math.pi/180)) 
    
    end = [end_lat,end_long]

    start = [start_lat,start_long]

    gmap_license="https://maps.googleapis.com/maps/api/js?key="+gmap_license

    Y_cord = (end[0]-start[0])/conversion
    X_cord = (end[1]-start[1])*(111*math.cos(start[0]*math.pi/180))

    if (Y_cord<=0 and X_cord<=0):
        theta =math.pi + math.atan(abs(Y_cord/X_cord))
    elif (Y_cord>=0 and X_cord<=0):
        theta =math.pi-math.atan(abs(Y_cord/X_cord))
    elif (Y_cord>=0 and X_cord>=0):
        theta =math.atan(abs(Y_cord/X_cord))
    elif (Y_cord<=0 and X_cord>=0):
        theta =2*math.pi-math.atan(abs(Y_cord/X_cord))


    
    # distance from start to end
    
    print(theta,Y_cord,X_cord)
    K = int(math.sqrt(X_cord**2 + Y_cord**2)/(l_scale*math.tan(alpha)))
    if K%2 ==0:
        K=K
    else:
        K=K+1
        print('layer is incremented by 1 to make it even')
        
    
    N= int(K/2)+1
    # edgelength
    L= math.sqrt((l_scale*math.tan(alpha))**2 + l_scale**2)
        # Get the grid in place
    Patchx=np.zeros((N,N))
    Patchy=np.zeros((N,N))
    Patchx,Patchy= Patchpoints(N, theta, alpha, l_scale)
    
    # Patchz holds the information of cost of constructon in the viscinity of vertices on the grid 
    Patchz=np.zeros((N,N))
    
    # 111 km per degree latitude change and 111.17km for per degree longitude change 
    #theta=-math.pi/3

    
    
    # Convert the nodes x,y location into latitude longitude information.
    
    k=0
    Lat=[1]*N**2
    Long=[1]*N**2
    for i in range(0,N):
        for j in range(0,N):
            Lat[k]=Patchy[i,j]*conversion+start[0]       
            Long[k]=Patchx[i,j]/(111.32*math.cos(Lat[k]*math.pi/180))+start[1]
            k=k+1

    
    # this is the path difference we can get at every drop or this is the minimum path 
    # separation we can achieve using this algorithm.

    delta=distance(Patchx[0,1],Patchy[0,1],Patchx[1,0],Patchy[1,0])
    d_half=int(d/(2*l_scale)) # layer to facilitate d separation
    buffer= int(request.POST['Buffer_layers']) # User can choose any number of layer that will be added to the 
                    #minimum separation constraint in order to make a choice of two 
                    #path combination out of many possible combinations.
    
    truncation_layer = int(d_half + buffer)

    truncation_layer_midpoint_start = [(Patchx[truncation_layer,0] + Patchx[0,truncation_layer])/2 , \
        (Patchy[truncation_layer,0] + Patchy[0,truncation_layer])/2]  

    # this is the distance of chord (truncation layer) from the circle.
    dr_truncated_layer = float(request.POST['Truncation_layerto_arc'])
    print(truncation_layer)
    # calculate the radius to restrict the last arc within dr_truncated_layer
    # the user does not need to enter it through GUI.          
    x=distance(truncation_layer_midpoint_start[0],\
                truncation_layer_midpoint_start[1],Patchx[truncation_layer,0],\
                    Patchy[truncation_layer,0])
    radius = (dr_truncated_layer**2 + x**2) / (2*dr_truncated_layer)      
    print(radius,x)
            
    
    # We define the center of the radial grid at the termnals.
    
   

    centery_start = truncation_layer_midpoint_start[1] + (radius - dr_truncated_layer) *math.sin(math.pi+theta)        
    centerx_start = truncation_layer_midpoint_start[0] + (radius - dr_truncated_layer) * math.cos(theta+(math.pi))       
    center_start_lat = start[0] + centery_start*conversion
    center_start_long = start[1] +centerx_start/(111.32*math.cos(center_start_lat*math.pi/180))

    # This decides number of rows in the radial grid 
    incremental_radius_start=float(request.POST['increment_start'])
    incremental_radius_end=float(request.POST['increment_end'])        
    
    # per delta angle 
    dtheta=math.asin(delta/(radius))  

    # number of nodes in radial grid 
    n_circum_start=int(2*math.pi/(dtheta) + 1)
    n_radial_start=int(radius/incremental_radius_start+1)

    # calculations for the end
    truncation_layer_midpoint_end = [(Patchx[N-1,N-1-truncation_layer] + Patchx[N-1-truncation_layer,N-1])/2 ,\
            (Patchy[N-1,N-1-truncation_layer] + Patchy[N-1-truncation_layer,N-1])/2]
    
    centery_end = truncation_layer_midpoint_end[1] + (radius - dr_truncated_layer) * math.sin(theta)       
    centerx_end= truncation_layer_midpoint_end[0] + (radius - dr_truncated_layer) * math.cos(theta)       
    center_end_lat = start[0] + centery_end * conversion        
    center_end_long= start[1] +centerx_end/(111.32*math.cos(center_end_lat*math.pi/180))
    
    n_circum_end=int(2*math.pi/(dtheta)+1)        
    n_radial_end=int(radius/incremental_radius_end+1)

# get the radial grid cordinates and latitude longitude.  

    radial_nodes_lat_start, radial_nodes_long_start,\
    radial_nodes_starty, radial_nodes_startx,radial_nodes_lat_end, radial_nodes_long_end,\
    radial_nodes_endy, radial_nodes_endx\
                = Patch_terminal(center_start_lat,center_start_long,center_end_lat,center_end_long,\
                    n_circum_start, n_radial_start,theta,radius, n_circum_end,\
                        n_radial_end,conversion)

    # convert it to list for sending it through post request to HTTP rendering.

    radial_nodes_lat_start=radial_nodes_lat_start.tolist()

    radial_nodes_long_start=radial_nodes_long_start.tolist()

    radial_nodes_starty=radial_nodes_starty.tolist()

    radial_nodes_startx=radial_nodes_startx.tolist()

    radial_nodes_lat_end=radial_nodes_lat_end.tolist()

    radial_nodes_long_end=radial_nodes_long_end.tolist()

    radial_nodes_endy=radial_nodes_endy.tolist()
   

    radial_nodes_endx=radial_nodes_endx.tolist()
   
    
    Patchx=Patchx.tolist()
    
    Patchy=Patchy.tolist()    
  
    request.session['dtheta'] = dtheta
    request.session['n_radial_start']=n_radial_start
    request.session['n_radial_end']=n_radial_end
    request.session['n_circum_start']=n_circum_start
    request.session['n_circum_end']=n_circum_end
    request.session['centerx_start'] = centerx_start
    request.session['centery_start'] = centery_start
    request.session['centerx_end'] = centerx_end
    request.session['centery_end'] = centery_end

    request.session['center_start_lat'] = center_start_lat
    request.session['center_start_long'] = center_start_long
    request.session['center_end_lat'] = center_end_lat
    request.session['center_end_long'] = center_end_long
    request.session['L'] = L
    request.session['theta'] = theta
  

    request.session['d'] = d
    request.session['incremental_radius_start'] = incremental_radius_start
    request.session['incremental_radius_end'] = incremental_radius_end        
    request.session['truncation_layer'] = truncation_layer
    request.session['N'] = N
    request.session['alpha'] = alpha
    request.session['gmap_license']=gmap_license
    request.session['K'] = K     
    request.session['start'] = start
    request.session['Patchx'] = Patchx
    request.session['Patchy'] = Patchy
    request.session['Lat'] = Lat
    request.session['Long'] = Long
    request.session['radius']=radius
    

    request.session['radial_nodes_long_start']=radial_nodes_long_start
    request.session['radial_nodes_lat_start']=radial_nodes_lat_start
    
    request.session['radial_nodes_starty']=radial_nodes_starty
    request.session['radial_nodes_startx']=radial_nodes_startx
    request.session['radial_nodes_lat_end']=radial_nodes_lat_end
    request.session['radial_nodes_long_end']=radial_nodes_long_end
    request.session['radial_nodes_endy']=radial_nodes_endy
    request.session['radial_nodes_endx']=radial_nodes_endx
    request.session['l_scale']= l_scale   

    return JsonResponse ({'radial_nodes_long_start' : radial_nodes_long_start,'truncation_layer':truncation_layer,\
            'radial_nodes_lat_start': radial_nodes_lat_start,'radial_nodes_long_end' : radial_nodes_long_end,\
            'radial_nodes_lat_end': radial_nodes_lat_end,'N':N,'K':K,'d':d,'start':start,'Patchx':Patchx,'Lat':Lat,\
        'Long':Long,'Patchy':Patchy,'gmap_license':gmap_license,'n_circum_start':n_circum_start,\
            'n_radial_start':n_radial_start,'n_circum_end':n_circum_end,'n_radial_end':n_radial_end})



def radial_SP(request):    
   
    dtheta = request.session['dtheta']
    n_radial_start=request.session['n_radial_start']
    n_radial_end=request.session['n_radial_end']
    n_circum_start=request.session['n_circum_start']
    n_circum_end=request.session['n_circum_end']
    truncation_layer = request.session['truncation_layer']

    radial_nodes_lat_start = request.session['radial_nodes_lat_start']
    radial_nodes_lat_start = np.array(radial_nodes_lat_start)

    radial_nodes_long_start=request.session['radial_nodes_long_start']
    radial_nodes_long_start=np.array(radial_nodes_long_start)

    radial_nodes_starty=request.session['radial_nodes_starty']
    radial_nodes_starty=np.array(radial_nodes_starty)

    radial_nodes_startx=request.session['radial_nodes_startx']
    radial_nodes_startx=np.array(radial_nodes_startx)

    
    radial_nodes_lat_end=request.session['radial_nodes_lat_end']
    radial_nodes_lat_end=np.array(radial_nodes_lat_end)

    radial_nodes_long_end=request.session['radial_nodes_long_end']
    radial_nodes_long_end=np.array(radial_nodes_long_end)

    radial_nodes_endy=request.session['radial_nodes_endy']
    radial_nodes_endy=np.array(radial_nodes_endy)

    radial_nodes_endx=request.session['radial_nodes_endx']
    radial_nodes_endx=np.array(radial_nodes_endx)
    l_scale = request.session['l_scale']
    alpha=request.session['alpha']
    radius = request.session['radius']
    theta = request.session['theta']
    d=request.session['d']
    K=request.session['K']
    N=request.session['N']
    start=request.session['start']
    incremental_radius_start= request.session['incremental_radius_start']
    incremental_radius_end= request.session['incremental_radius_end']
    center_start_lat = request.session['center_start_lat'] 
    center_start_long = request.session['center_start_long'] 
    center_end_lat =  request.session['center_end_lat']
    center_end_long = request.session['center_end_long']
  
    Longitude_start=[]
    Latitude_start=[]
    Cost_start=[]

    Longitude_end=[]
    Latitude_end=[]
    Cost_end=[]

    # user can define a start point anywhere within the radial grid.
   
    # start_lat = float(request.POST['start_lat'])
    # start_long = float(request.POST['start_long'])
    # end_lat = float(request.POST['end_lat'])
    # end_long = float(request.POST['end_long'])

    # # we check if the user has given the lat-long within the radial graph or not
    # conversion=1/111

    # centr_to_starty = (start_lat - center_start_lat)/conversion
    # centr_to_startx = (start_long - center_start_long) * \
    #     (111.32*math.cos(center_start_lat*math.pi/180))

    # angle_start = math.atan(centr_to_starty/centr_to_startx)
    
    
    # if math.sqrt(centr_to_starty**2 + centr_to_startx**2) > radius:
        
    #     start_lat = center_start_lat + radius* math.cos(angle_start-theta)*conversion
    #     start_long = center_start_long + radius* math.sin(angle_start-theta)/(111.32*math.cos(center_start_lat*math.pi/180))
    #     centr_to_starty = (start_lat - center_start_lat)/conversion
    #     centr_to_startx = (start_long - center_start_long) * (111.32*math.cos(center_start_lat*math.pi/180))
    #     print('The start is shifted within the radial grid to', start_lat,start_long)

    # centr_to_endy = (end_lat - center_end_lat)/conversion
    # centr_to_endx = (end_long - center_end_long) * \
    #     (111.32*math.cos(center_end_lat*math.pi/180))

    # angle_end = math.atan(centr_to_endy/centr_to_endx)

    # if math.sqrt(centr_to_endy**2 + centr_to_endx**2) > radius:        
    #     end_lat = center_end_lat + radius* math.cos(angle_end+theta)*conversion
    #     end_long = center_end_long + radius* math.sin(angle_end+theta)/(111.32*math.cos(center_end_lat*math.pi/180))
    #     centr_to_endy = (end_lat - center_end_lat)/conversion
    #     centr_to_endx = (end_long - center_end_long) * \
    #         (111.32*math.cos(center_end_lat*math.pi/180))
        
    #     print('The end is shifted within the radial grid to', end_lat,end_long)

    
    # compute start node number on the radial grid
    start_node_number=0
    #int(angle_start/dtheta) +\
    #       n_circum_start*int(math.sqrt(centr_to_starty**2 + centr_to_startx**2)/incremental_radius_start)

    end_node_number = 0
    # int(angle_end/dtheta) +\
    #         n_circum_end*int(math.sqrt(centr_to_endy**2 + centr_to_endx**2)/incremental_radius_end)

    # Compute the target node number on the radial grid.
    last_angle_start=math.asin(truncation_layer*l_scale/(radius))

    last_column_left=int(last_angle_start/dtheta)
    print(last_column_left)
    target_nodes_start_j1 = list(np.linspace(last_column_left,0,last_column_left+1))
    
    target_nodes_start_j = target_nodes_start_j1+[(n_circum_start-1-i) for i in np.linspace(1,last_column_left,last_column_left)]
    target_nodes_start_i = n_radial_start-1
    target_nodes_start = [(target_nodes_start_i) * n_circum_start + target_nodes_start_j[i] for i in range(len(target_nodes_start_j))]
    target_nodes_start = target_nodes_start[::-1] # for matching the indexing convention of
    print(target_nodes_start)   
        
    #target_nodes_start=[N_start-1,N_start-2]
    last_angle_end=math.asin(l_scale*truncation_layer/(radius))
    last_column_left=int(last_angle_end/dtheta)
    print(last_column_left)
    target_nodes_end_j1 = list(np.linspace(last_column_left,0,last_column_left+1))
    target_nodes_end_j = target_nodes_end_j1 + [(n_circum_start-1-i) for i in np.linspace(1,last_column_left,last_column_left)]
    target_nodes_end_i = n_radial_end - 1
    target_nodes_end = [(target_nodes_end_i) * n_circum_end + target_nodes_end_j[i] for i in range(len(target_nodes_end_j))]
    
    # if  request.is_ajax() and request.POST:
    nTimes=int(request.POST['Grid_density'])
    myfile = request.FILES["radial_start"]
    # reader = csv.DictReader(str(myfile.name))
    # print(os.path.getsize(reader))
    # #lines=myfile.read()
    # #print(len(lines))
    # #Radial_Start=myfile.readlines()
    # for row in reader:       
    #     Latitude_start.append(float(row['Latitude']))
    #     print(Latitude_start)
    #     Longitude_start.append(float(row['Longitude']))
    #     Cost_start.append(float(row['Cost']))

    # print(Longitude_start)


    fs = FileSystemStorage() #defaults to   MEDIA_ROOT  
    ABC = fs.save('radial_start.csv', myfile)

    with open('media/radial_start.csv', 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            Latitude_start.append(float(row['Latitude']))
            Longitude_start.append(float(row['Longitude']))
            Cost_start.append(float(row['Cost']))
  

   
    myfile1 = request.FILES["radial_end"]
    fs = FileSystemStorage()
    CFD = fs.save('radial_end.csv', myfile1) 
    
  
    with open('media/radial_end.csv','r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            Latitude_end.append(float(row['Latitude']))
            Longitude_end.append(float(row['Longitude']))
            Cost_end.append(float(row['Cost']))
    

    
    #Interpolate costs at nodes
    # Inverse distance weight algorithm. We get a list of interpolated points for the unknown xi, yi,
    # using known points x, y, z.
    # 
    Lat_lim = float(request.POST['Lat_limit'])
    Long_lim = float(request.POST['Long_limit'])
    
    Cost_start_input = IDW(np.array(Latitude_start),np.array(Longitude_start),np.array(Cost_start),radial_nodes_lat_start.\
        reshape(int(n_circum_start)*n_radial_start,), radial_nodes_long_start.reshape((n_circum_start)*n_radial_start,),Lat_lim,Long_lim)

    

    Cost_end_input=IDW(np.array(Latitude_end),np.array(Longitude_end),np.array(Cost_end),radial_nodes_lat_end.\
        reshape(int(n_circum_end)*n_radial_end,), radial_nodes_long_end.reshape((n_circum_end)*n_radial_end,),Lat_lim,Long_lim)

    
    Patchz_start=np.array(Cost_start_input).reshape(n_radial_start,n_circum_start)
    print(Patchz_start[1,:])

    Patchz_end=np.array(Cost_end_input).reshape(n_radial_end,n_circum_end)
    # Check for the inputs to be within 
    if min(Cost_start_input)==math.inf:
        print('There is some mistake in data input')
    if min(Cost_end_input)==math.inf:
        print('There is some mistake in data input')
                    

    r_start_circum=(nTimes-1)*(n_circum_start-1)+1
    r_start_radial=(nTimes-1)*(n_radial_start-1)+1

    r_end_circum=(nTimes-1)*(n_circum_end-1)+1
    r_end_radial=(nTimes-1)*(n_radial_end-1)+1

    zvals_start=np.zeros((r_start_radial,r_start_circum))
    zvals_end=np.zeros((r_end_radial,r_end_circum))    

    zvals_start = interpolation_radial(nTimes,zvals_start,Patchz_start)
    zvals_end = interpolation_radial(nTimes,zvals_end,Patchz_end) 
    

    # edges costs of radial graphs,

    Cost_radial_start = Cost_radial(n_circum_start,zvals_start,nTimes,n_radial_start,\
                                radial_nodes_startx, radial_nodes_starty)

    Cost_radial_end = Cost_radial(n_circum_end,zvals_end,nTimes,n_radial_end,\
                                radial_nodes_endx, radial_nodes_endy)

    #print(Cost_radial_start[target_nodes_start[0]//n_circum_start,target_nodes_start[0]%n_circum_start],  Cost_radial_end[target_nodes_end[0]//n_circum_end,target_nodes_end[0]%n_circum_end])


    # total number of nodes in the radial grid.

    N_start = int(n_circum_start*n_radial_start)
    N_end  = int(n_circum_end*n_radial_end)

    # start computing the radial SP
    highvalue=int(request.POST['High_val'])
    sPath_terminal_start=[]
    sPath_terminal_end=[]

    
    # check if the targets are not sitting at infinity cost.
    for i in range(len(target_nodes_start)):
        res = all(ele == math.inf for ele in Cost_radial_start[int(target_nodes_start[i]),:])
        
        if res==True:
            print('shift the target, its in infinite cost zone.')


    for i in range(len(target_nodes_start)):

        results = bellman_ford(Cost_radial_start,N_start,start_node_number,int(target_nodes_start[i]),highvalue)
        
        sPath_terminal_start.append(results)
    
       

    for i in range(len(target_nodes_end)):        
        sPath_terminal_end.append(bellman_ford(Cost_radial_end,N_end,end_node_number,target_nodes_end[i],highvalue))

    print(sPath_terminal_start,sPath_terminal_end)
    # get list of costs from source to targets in truncation layer.
    truncated_layer_start_cost = [sPath_terminal_start[i][1] for i in range(len(sPath_terminal_start))]
    truncated_layer_end_cost = [sPath_terminal_end[i][1] for i in range(len(sPath_terminal_end))]
    
   
    
    radial_nodes_long_end=radial_nodes_long_end.tolist()
    radial_nodes_long_start=radial_nodes_long_start.tolist()
    radial_nodes_lat_end=radial_nodes_lat_end.tolist()
    radial_nodes_lat_start=radial_nodes_lat_start.tolist()

    request.session['truncated_layer_start_cost'] = truncated_layer_start_cost
    request.session['truncated_layer_end_cost'] = truncated_layer_end_cost

    with open('sPath_terminal.txt','w') as file_list:
        for lines in sPath_terminal_start:
            file_list.write('%s\n' % lines)
    with open('sPath_terminal_end.txt','w') as file_list:
        for lines in sPath_terminal_end:
            file_list.write('%s\n' % lines)

    os.remove(os.path.join(settings.MEDIA_ROOT, ABC))
    os.remove(os.path.join(settings.MEDIA_ROOT, CFD))
    

    return JsonResponse ({'sPath_terminal_start':sPath_terminal_start,'sPath_terminal_end':sPath_terminal_end})

def shortest_path1(request):
    Lat=request.session['Lat']
    Long=request.session['Long']
    alpha=request.session['alpha']
    N=request.session['N']
    d=request.session['d']
    K=request.session['K']
    Lat=request.session['Lat']
    L = request.session['L']
    Long=request.session['Long']
    Patchx=request.session['Patchx']
    Patchy=request.session['Patchy']   
    start=request.session['start']
    truncated_layer_start_cost= np.array(request.session['truncated_layer_start_cost'])
    truncated_layer_end_cost= np.array(request.session['truncated_layer_end_cost'])
    truncation_layer = request.session['truncation_layer']
    center_start_lat = request.session['center_start_lat'] 
    center_start_long = request.session['center_start_long'] 
    center_end_lat =  request.session['center_end_lat']
    center_end_long = request.session['center_end_long']
    l_scale = request.session['l_scale']

    n=int(1/48*((K+2)*(K+4)*(2*K+6)))
    Patchz=np.zeros((N,N))   
    Patchx=np.array(Patchx)
    Patchy=np.array(Patchy)
# surface patch density and edgelength
   
    nTimes=int(request.POST['Grid_density'])
    myfile = request.FILES["Cost_matrix"]
    fs = FileSystemStorage()
    Cost_matrix = fs.save('Cost_matrix.csv', myfile)    

    L_Diag = 2*L*math.sin(alpha)  
    Longitude_input = []
    Latitude_input = []
    Cost_input = []  
    

    # Reading the file from filestorage
    with open('media/Cost_matrix.csv', 'w', newline='') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            Latitude_input.append(float(row['Latitude']))
            Longitude_input.append(float(row['Longitude']))
            Cost_input.append(float(row['Cost']))

    

   
    Latitude_input=np.array(Latitude_input)
    Longitude_input=np.array(Longitude_input)
    Cost_input=np.array(Cost_input)
    Lat=np.array(Lat)
    Long=np.array(Long)
    
    #Interpolate costs at nodes
    # Inverse distance weight algorithm. We get a list of interpolated points for the unknown xi, yi,
    # using known points x, y, z.

    Lat_lim = float(request.POST['Lat_limit'])
    Long_lim = float(request.POST['Long_limit'])

    Cost=IDW(Latitude_input,Longitude_input,Cost_input,Lat,Long,Lat_lim, Long_lim)

    Hash=np.zeros((n,3))
    Hash=hash1(n,K)

    N_K=np.array([1.0]*(K+1)) # number of nodes in each layer Note that there are K layers
                                                            #between K+1 planes.
    i=1
    
    # Number of nodes in each layer. Each layer has has number of nodes as sum of integers from zero to that layer.
    
    while i<=K+1:
        if i<=K/2:
            N_K[i-1]=int(i*(i+1)/2)
        else:
            N_K[i-1]=int((K-i+2)*(K-i+3)/2)
        i=i+1

    # Writing the interpolated node to the media folder for the user to view manipulate and update the cost file if required.
     
    with open('media/diamond_frame_nodes.csv', "w") as csv_file:
        writer = csv.writer(csv_file, delimiter=',',lineterminator='\n')
        writer.writerow(['Latitude','Longitude','Cost'])
        for i in range(0, len(Lat)):
            writer.writerow([Lat[i],Long[i],Cost[i]])

    k=0    
    for i in range(0,N):
        for j in range(0,N):
            Patchz[i,j]=Cost[k]
            k=k+1        
    r=(nTimes-1)*(N-1)+1
    
    zvals=np.zeros((r,r))
    zvals=interpolation(nTimes,N,Patchz)

    Cost_horizontal,Cost_vertical = Cost_diamondgraph(N,zvals,nTimes,L,L_Diag)
    
    #d_scaled = int(d/l_sc

    Node_val=node_val(Hash,Cost_horizontal,Cost_vertical,N_K,truncated_layer_start_cost, truncated_layer_end_cost,K,d,truncation_layer)

  
    # Shortest path computation           
    sPath=np.zeros((K,3))    
    sPath[K-1,:]=Node_val[n-1,:]
    for m in range(K-2,-1,-1):
        sPath[m,:]=Node_val[int(sPath[m+1,1]),:]
    print(Cost_horizontal[0,:],sPath)

    request.session['N_K']=N_K.tolist()
    
    sPath=sPath.tolist()
    Hash=Hash.tolist()
    request.session['Hash'] = Hash
    request.session['N'] = N
    

    Cost_horizontal = np.where(Cost_horizontal == math.inf , 1000000, Cost_horizontal)
    Cost_vertical = np.where(Cost_vertical == math.inf , 1000000, Cost_vertical)
    

    Cost_horizontal = Cost_horizontal.tolist()
    Cost_vertical = Cost_vertical.tolist()
  
   
    
    Patchx = Patchx.tolist()
    Patchy = Patchy.tolist()
    Lat = Lat.tolist()
    Long = Long.tolist()    

    request.session['Cost_horizontal'] = Cost_horizontal
    request.session['Cost_vertical'] = Cost_vertical
    os.remove(os.path.join(settings.MEDIA_ROOT, Cost_matrix))
    
    return JsonResponse({'sPath':sPath,'start': start,'Hash':Hash,'Patchx':Patchx,'Patchy':Patchy,'Lat':Lat,\
        'Long':Long,'Cost_horizontal':Cost_horizontal,'Cost_vertical':Cost_vertical,'K':K,'d':d, 'N':N,\
             'truncation_layer':truncation_layer, 'center_start_lat':center_start_lat, \
                'center_start_long':center_start_long,'center_end_lat':center_end_lat, 'center_end_long':center_end_long})

def Modify_cost(request):
# request is a REST api framework using which we can request variables saved as session in other views modules...
    Cost_horizontal=np.array(request.session['Cost_horizontal'])
    Cost_vertical=np.array(request.session['Cost_vertical'])
    
    K=request.session['K']
    
    d=request.session['d']
    truncated_layer_start_cost = np.array(request.session['truncated_layer_start_cost'])
    truncated_layer_end_cost = np.array(request.session['truncated_layer_end_cost'])
    truncation_layer =request.session['truncation_layer']
    
    Patchx=request.session['Patchx']
    Patchy=request.session['Patchy']
    N_K=request.session['N_K']
    N_K=np.array(N_K)  
    N=request.session['N']    
    K=request.session['K']
    d=request.session['d']
    start=request.session['start']
    center_start_lat = request.session['center_start_lat'] 
    center_start_long = request.session['center_start_long'] 
    center_end_lat =  request.session['center_end_lat']
    center_end_long = request.session['center_end_long']

    Hash = request.session['Hash']     
    Hash=np.array(Hash) 
    
    n=int(1/48*((K+2)*(K+4)*(2*K+6)))

    # receive list of nodes whose value is infinity
    
    node_ids=  request.POST.getlist('nodeid_forcostchange[]')
    COST =  request.POST.getlist('COST[]')
    
    if len(node_ids)!=0:
        for nodes,costs in zip(node_ids,COST):            
            i=int(nodes)
            j=costs.split(',')
            print(costs)            
            Cost_horizontal[i//N,i%N +1] = float(j[0])            
            Cost_vertical[i//N,i%N + 1] = float(j[1])
            

    
    sPath1=np.zeros((K,3))
    
    Node_val= node_val(Hash,Cost_horizontal,Cost_vertical,N_K,truncated_layer_start_cost, truncated_layer_end_cost,K,d,truncation_layer)
    sPath1[K-1,:]=Node_val[n-1,:]
    
    for m in range(K-2,-1,-1):
        sPath1[m,:]=Node_val[int(sPath1[m+1,1]),:]
 
    sPath1=sPath1.tolist()
    Hash=Hash.tolist()    
    N_K=N_K.tolist()
    Cost_horizontal=Cost_horizontal.tolist()
    Cost_vertical=Cost_vertical.tolist()
    
   

    return JsonResponse({'Cost_horizontal':Cost_horizontal,'Cost_vertical':Cost_vertical,'sPath1':sPath1,'start': start,'K':K,'d':d,'Hash':Hash,'Patchx':Patchx,'Patchy':Patchy,'center_start_lat':center_start_lat,'center_start_long':center_start_long,'center_end_lat':center_end_lat, 'center_end_long':center_end_long})

    

    
     
    



